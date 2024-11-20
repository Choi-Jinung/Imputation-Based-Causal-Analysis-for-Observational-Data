library(dplyr)

# Implementation of variance estimation introduced in "Statistical matching using fractional imputation (Jae Kwang Kim el al., 2016)"

# Argument
## 1. FI_single_output: output of the function 'FI_single'
variance_estimation_single <- function(FI_single_output) {
  gamma <- FI_single_output$gamma; kappa <- gamma[-length(gamma)]; sig2 <- gamma[length(gamma)] 
  phi <- FI_single_output$est_para; beta <- phi[-length(phi)]; tau2 <- phi[length(phi)]; 
  y1_bar_A <- FI_single_output$observed_bar_prop; y1_bar_hat_B <- FI_single_output$observed_bar; y1_bar_hat <- FI_single_output$observed_bar_total 
  dataA <- FI_single_output$prop_data; dataB <- FI_single_output$target_data; imputed_dataB <- FI_single_output$imputed_data
  n_A <- nrow(dataA); n_B <- nrow(dataB)
  proposal_fit <- FI_single_output$proposal_fit; IV <- FI_single_output$IV
  M_B <- FI_single_output$M_B
  s2_predictors <- FI_single_output$z_names
  s1_predictors <- colnames(proposal_fit$model)[-1]
  p_s1 <- length(s1_predictors) + 2
  y1_names <- FI_single_output$y1_names; y2_names <- FI_single_output$y2_names
  
  # Building Blocks
  XX_s2 <- as.matrix(cbind(1, imputed_dataB[, s2_predictors]))
  XX_s1 <- as.matrix(cbind(1, imputed_dataB[, s1_predictors]))
  
  scores2 <- cbind(sweep(XX_s2, MARGIN = 1, STATS = (tau2^-1)*(imputed_dataB[, y2_names] - XX_s2%*%beta), FUN = "*"), -(tau2^-1)/2 + ((tau2^-2)/2)*((imputed_dataB[, y2_names] - XX_s2%*%beta)^2))
  scores1 <- cbind(sweep(XX_s1, MARGIN = 1, STATS = (sig2^-1)*(imputed_dataB[, y1_names] - XX_s1%*%kappa), FUN = "*"), -(sig2^-1)/2 + ((sig2^-2)/2)*((imputed_dataB[, y1_names] - XX_s1%*%kappa)^2))
  
  scores2_bar <- data.frame(scores2, ID = imputed_dataB$id_B, frac_w = imputed_dataB$w)
  scores2_bar <- scores2_bar %>% group_by(ID) %>% summarise(across(everything(), ~weighted.mean(., frac_w))); scores2_bar$ID <- NULL; scores2_bar$frac_w <- NULL; scores2_bar <- as.matrix(scores2_bar)
  scores1_bar <- data.frame(scores1, ID = imputed_dataB$id_B, frac_w = imputed_dataB$w)
  scores1_bar <- scores1_bar %>% group_by(ID) %>% summarise(across(everything(), ~weighted.mean(., frac_w))); scores1_bar$ID <- NULL; scores1_bar$frac_w <- NULL; scores1_bar <- as.matrix(scores1_bar)
  
  B21 <- t(sweep(scores2, MARGIN = 1, STATS = imputed_dataB$w, FUN = "*"))%*%(scores1 - scores1_bar[imputed_dataB$id_B, ])
  B22 <- t(sweep(scores2, MARGIN = 1, STATS = imputed_dataB$w, FUN = "*"))%*%(scores2 - scores2_bar[imputed_dataB$id_B, ])
  
  I22_11 <- -(tau2^-1)*t(sweep(XX_s2, MARGIN = 1, STATS = sqrt(imputed_dataB$w), FUN = "*"))%*%sweep(XX_s2, MARGIN = 1, STATS = sqrt(imputed_dataB$w), FUN = "*")
  I22_12 <- -(tau2^-2)*t(sweep(XX_s2, MARGIN = 1, STATS = sqrt(imputed_dataB$w), FUN = "*"))%*%sweep((imputed_dataB[, y2_names] - XX_s2%*%beta), MARGIN = 1, STATS = sqrt(imputed_dataB$w), FUN = "*")
  I22_22 <- (n_B/2)*(tau2^-2) - (tau2^-3)*sum(imputed_dataB$w*(imputed_dataB[, y2_names] - XX_s2%*%beta)^2)
  I22_former <- rbind(cbind(I22_11, I22_12), cbind(t(I22_12), I22_22))
  I22 <- I22_former + B22
  
  V_S2_bar <- t(scores2_bar)%*%scores2_bar
  
  scores1_A <- cbind(sweep(cbind(1, proposal_fit$model[, -1]), MARGIN = 1, STATS = (sig2^-1)*proposal_fit$residuals, FUN = "*"), -(sig2^-1)/2 + ((sig2^-2)/2)*(proposal_fit$residuals^2))
  scores1_A <- as.matrix(scores1_A)
  V_S1 <- t(scores1_A)%*%scores1_A
  I11_11 <- sig2*t(as.matrix(cbind(1, proposal_fit$model[, -1])))%*%as.matrix(cbind(1, proposal_fit$model[, -1]))
  I11_12 <- matrix(apply(sweep(cbind(1, proposal_fit$model[, -1]), STATS = proposal_fit$residuals, MARGIN = 1, FUN = "*") , MARGIN = 2, FUN = sum), nrow = p_s1-1)
  I11_22 <- -n_A/2 + (sig2^-1)*sum(proposal_fit$residuals^2)
  I11 <- ((sig2^-2))*rbind(cbind(I11_11, I11_12), cbind(t(I11_12), I11_22))
  KK <- B21%*%solve(I11)
  V2 <- KK%*%V_S1%*%t(KK)
  
  V_beta <- solve(I22)%*%(V_S2_bar+V2)%*%t(solve(I22))
  lower.bound <- phi - sqrt(diag(V_beta))*qnorm(0.975)
  upper.bound <- phi + sqrt(diag(V_beta))*qnorm(0.975)
  CI_95 <- cbind(lower.bound, upper.bound)
  
  var_est_output <- list(V_beta=V_beta, CI_95=CI_95)
  return(var_est_output)
}
