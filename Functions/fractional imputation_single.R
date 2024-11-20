library(dplyr)
library(nleqslv)

# Implementation of "Statistical matching using fractional imputation (Jae Kwang Kim et al., 2016)".

# Arguments
## 1. dataA: data frame with X and "y1"
## 2. dataB: data frame with X and "y2". Predicting y1 values of this data set is the goal.
## 3. IV: string vector for names of non-response instrumental variables. 
## It specifies the names of variables among X to be excluded from f(y2 | X, y1).
## 4. M_B: imputation size
## 5. para_curr: initial guess of the parameter
## 6. EM: boolean varialbe indicating whether to use EM algorithm.
## If EM==F, the mean score equation is directly solved using the 'nleqslv' package. 

FI_single <- function(dataA, dataB, IV, M_B, para_curr, EM) {
  n_A <- nrow(dataA); n_B <- nrow(dataB)
  y1_names <- setdiff(colnames(dataA), colnames(dataB))
  y2_names <- setdiff(colnames(dataB), colnames(dataA))
  x_names <- intersect(colnames(dataA), colnames(dataB))
  z_names <- setdiff(c(y1_names, x_names), IV)
  
  # Initialization
  beta1_curr <- para_curr[-length(para_curr)] 
  var1_curr <- para_curr[length(para_curr)]
  
  # Step1: Generate M imputed values for an i-th unit in dataB
  
  # Fit the proposal distribution with dataA
  proposal_fit <- lm(as.formula(paste(y1_names, "~.")), data=dataA) 
  proposal_mean <- predict(proposal_fit, newdata = dataB)
  proposal_sd <- sqrt(sum(proposal_fit$residuals^2)/proposal_fit$df.residual)
  gamma <- c(proposal_fit$coefficients, proposal_sd^2)
  
  # Generate imputed values from the fitted proposal distribution
  # and augment the dataB with the imputed values
  imputed_values <- proposal_mean[rep(1:nrow(dataB), each = M_B)] + rnorm(n = nrow(dataB)*M_B, mean = 0, sd = proposal_sd)
  id_B <- rep(1:nrow(dataB), each = M_B)
  imputed_dataB <- cbind(id_B, dataB[rep(1:nrow(dataB), each = M_B), ], imputed_values)
  head(dataB)
  head(imputed_dataB, 50)
  colnames(imputed_dataB)[length(colnames(imputed_dataB))] <- y1_names
  
  # Use EM algorithm to solve the mean score equation
  if (EM) {
    tt <- 0
    X <- as.matrix(cbind(1, imputed_dataB[, c(z_names)])) # covariate
    Y <- imputed_dataB[, y2_names] # imputed values
    max_iter <- 1000
    conv.path <- matrix(NA, nrow=2*max_iter, ncol=length(z_names)+3)
    repeat {
      tt <- tt+1
      
      # Step2: Assign fractional weight to each imputed value (E-step)
      # Compute fractional weight
      imputed_dataB$w <- dnorm(x=imputed_dataB[, y2_names], mean=as.matrix(cbind(1, imputed_dataB[, z_names]))%*%beta1_curr, sd=sqrt(var1_curr)) 
      imputed_dataB$w <- imputed_dataB %>% group_by(id_B) %>% reframe(w = w/sum(w)) %>% .[[2]]
      
      # Step3: Solve the fractionally imputed score equation (M-step)
      W <- imputed_dataB[, "w"]
      fit_next <- lm(Y~-1+X, weights = W)
      para_next <- c(fit_next$coefficients, sum(W*fit_next$residuals^2)/sum(W))
      
      # Step4: Check convergence
      dist <- sqrt(sum((para_next-para_curr)^2))
      conv.path[tt, ] <- c(para_next, dist)
      # If theta_(t) and theta_(t+1) are close enough, break the iteration
      if (dist < 1e-6 | tt > max_iter) {
        cat("Iteration:", tt, "Distance:", dist, "\n") 
        break
      } # If else, update the current parameters
      else {
        cat("Iteration:", tt, "Distance:", dist, "\n")
        para_curr <- para_next
        beta1_curr <- para_curr[-length(para_curr)]
        var1_curr <- para_curr[length(para_curr)]
      }
    }
    conv.path <- conv.path[is.na(conv.path[, 1])==F, ]
    termcd <- NULL
    
    # If EM==F, the mean score equation will be directly solved rather than using EM algorithm
  } 
  # else {
  #   # Define the mean score function
  #   score2 <- function(phi) {
  #     beta1 <- phi[-length(phi)]; var1 <- phi[length(phi)] 
  #     imputed_dataB$w <- dnorm(x=imputed_dataB[, y2_names], mean=as.matrix(cbind(1, imputed_dataB[, z_names]))%*%beta1, sd=sqrt(var1)) 
  #     imputed_dataB$w <- imputed_dataB %>% group_by(id_B) %>% reframe(w = w/sum(w)) %>% .[[2]]
  #     mu_beta <- as.matrix(cbind(1, imputed_dataB[, z_names]))
  #     resid <- imputed_dataB[, y2_names] - mu_beta%*%beta1
  #     unit_scores <- cbind(sweep(mu_beta, MARGIN = 1, STATS = (var1^-1)*resid, FUN = "*"), -(var1^-1)/2 + ((var1^-2)/2)*(resid^2))
  #     weighted_unit_scores <- sweep(unit_scores, MARGIN = 1, STATS = imputed_dataB$w, FUN = "*")
  #     score2_out <- apply(weighted_unit_scores, MARGIN = 2, FUN = sum)/(n_A+n_B)
  #     return(score2_out)
  #   }
  #   
  #   # Solve the mean score equation
  #   nleqslv_result <- nleqslv(para_curr, score2); termcd <- nleqslv_result$termcd
  #   para_next <- nleqslv_result$x
  #   beta1 <- para_next[-length(para_next)]
  #   var1 <- para_next[length(para_next)]
  #   
  #   # Calculate the fractional weights for each imputed value
  #   imputed_dataB$w <- dnorm(x=imputed_dataB[, y2_names], mean=as.matrix(cbind(1, imputed_dataB[, z_names]))%*%beta1, sd=sqrt(var1)) 
  #   imputed_dataB$w <- imputed_dataB %>% group_by(id_B) %>% reframe(w = w/sum(w)) %>% .[[2]]
  #   conv.path <- NULL; fit_next <- NULL
  # }
  
  imputed_values_id <- cbind(imputed_values, imputed_dataB[, c("id_B", "w")])
  missing_bar <- weighted.mean(imputed_dataB[, y1_names], imputed_dataB$w) # mean of y1 in dataB
  observed_bar <- weighted.mean(imputed_dataB[, y2_names], imputed_dataB$w) # mean of y2 in dataB
  
  result <- list(gamma = gamma, est_para=para_next, 
                 missing_bar = missing_bar, observed_bar = observed_bar,
                 prop_data = dataA, target_data = dataB, imputed_data = imputed_dataB, 
                 proposal_fit = proposal_fit, fit_next = fit_next, IV = IV, M_B = M_B, conv.path=conv.path,
                 z_names = z_names, termcd = termcd, 
                 y1_names = y1_names, y2_names = y2_names)
                 
  return(result)
}
