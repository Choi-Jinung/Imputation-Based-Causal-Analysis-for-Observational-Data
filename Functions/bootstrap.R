library('doParallel')
library('doRNG')

source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/fractional imputation_single.R")
source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/Regression_CI.R")
source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/PSM.R")

# source("/Users/a82103/Desktop/Research/Causal Inference/Code_final/Functions/fractional imputation_single.R")
# source("/Users/a82103/Desktop/Research/Causal Inference/Code_final/Functions/Regression_CI.R")
# source("/Users/a82103/Desktop/Research/Causal Inference/Code_final/Functions/PSM.R")

# Bootstrap variance estimation of ATT, ATC and ATE

# Arguments
## 1. control: control group data
## 2. treatment: treatment group data
## 3. IV: string vector for names of non-response instrumental variables. 
## It specifies the names of variables among X to be excluded from f(y2 | X, y1).
## 4. M_B: imputation size
## 5. para_curr: initial guess of the parameter
## 6. EM: boolean varialbe indicating whether to use EM algorithm.
## If EM==F, the mean score equation is directly solved using the 'nleqslv' package. 
## 7. covariate: predictor variables to be used in the CIA regression
## 8. B: bootstrap sample size
## 9. seed: random seed for reproduction

var_est<- function(control, treatment, IV, M_B, para_curr, EM, covariate, B, seed) {
  start <- proc.time()
  
  n1 <- nrow(treatment); n0 <- nrow(control)
  y1_names <- setdiff(colnames(treatment), colnames(control))
  y0_names <- setdiff(colnames(control), colnames(treatment))
  control[, y1_names] <- NA
  treatment[, y0_names] <- NA
  merged_data <- rbind(control, treatment)
  control[, y1_names] <- NULL
  treatment[, y0_names] <- NULL
  
  # Obtain CIA, PSM and FI estimates for one bootstrap sample
  boot_one_time <- function(merged_data, y1_names, y0_names, IV, M_B, para_curr, EM, covariate, it) {
    
    # Generate a bootstrap sample
    boot_sampling <- function(data) {
      n <- nrow(data)
      index <- sample(1:n, n, replace=T)
      resampled_data <- data[index, ]
      return(resampled_data)
    }
    merged_data_boot <- boot_sampling(merged_data)
    data_psm_boot <- merged_data_boot %>% mutate(t = case_when(is.na(y0)~1, is.na(y1)~0))
    data_psm_boot <- data_psm_boot %>% mutate(y = case_when(is.na(y0)~y1, is.na(y1)~y0))
    data_psm_boot$y0 <- NULL; data_psm_boot$y1 <- NULL
    
    # Split the bootstrap sample into control group and treatment group
    split_data_boot <- function(boot_data, y1_names, y0_names) {
      control_boot <- boot_data[is.na(boot_data[, y1_names]), ]; control_boot$y1 <- NULL
      treatment_boot <- boot_data[is.na(boot_data[, y0_names]), ]; treatment_boot$y0 <- NULL
      split_boot_data_out <- list(control_boot=control_boot, treatment_boot=treatment_boot)
      return(split_boot_data_out)
    }
    split_data_boot_list <- split_data_boot(merged_data_boot, y1_names, y0_names)
    control_boot <- split_data_boot_list$control_boot; treatment_boot <- split_data_boot_list$treatment_boot
    n0_boot <- nrow(control_boot); n1_boot <- nrow(treatment_boot)
    
    # Propensity Score Matching (PSM)
    PSM_ATT <- PSM(data=data_psm_boot, 
                   method="full", estimand = "ATT", estimand_group = 1,
                   covariates_ps = c("x1","x2","x3","x4"),
                   covariates_out = c("x1","x2","x3","x4"))
    PSM_ATC <- PSM(data=data_psm_boot, 
                   method="full", estimand = "ATC", estimand_group = 0,
                   covariates_ps = c("x1","x2","x3","x4"),
                   covariates_out = c("x1","x2","x3","x4"))
    PSM_ATE <- PSM(data=data_psm_boot, 
                   method="full", estimand = "ATE", estimand_group = "total",
                   covariates_ps = c("x1","x2","x3","x4"),
                   covariates_out = c("x1","x2","x3","x4"))
    PSM_ATE_estimates <- c(PSM_ATT$est$estimate, PSM_ATC$est$estimate, PSM_ATE$est$estimate)
    
    # Regression Imputation (CIA)
    reg_ATT_list <- Regression_CI(control_boot, treatment_boot, covariate)
    reg_ATC_list <- Regression_CI(treatment_boot, control_boot, covariate)
    reg_ATT <- reg_ATT_list$observed_bar - reg_ATT_list$missing_bar
    reg_ATC <- reg_ATC_list$missing_bar - reg_ATC_list$observed_bar
    reg_ATE <- (n1_boot*reg_ATT + n0_boot*reg_ATC)/(n1_boot+n0_boot)
    reg_ATE_estimates <- c(reg_ATT, reg_ATC, reg_ATE)
    
    # Fractional Imputation
    impute_boot_control <- FI_single(dataA=treatment_boot, dataB=control_boot, IV=IV, M_B=M_B, para_curr=para_curr, EM=T)
    impute_boot_treatment <- FI_single(dataA=control_boot, dataB=treatment_boot, IV=IV, M_B=M_B, para_curr=para_curr, EM=T)
    FI_ATC <- impute_boot_control$missing_bar - impute_boot_control$observed_bar
    FI_ATT <- impute_boot_treatment$observed_bar - impute_boot_treatment$missing_bar
    FI_ATE <- (n0_boot*FI_ATC + n1_boot*FI_ATT)/(n0_boot+n1_boot)
    FI_ATE_estimates <- c(FI_ATT, FI_ATC, FI_ATE)
    
    boot_one_time_output <- list(PSM_ATE_estimates=PSM_ATE_estimates,
                                 reg_ATE_estimates=reg_ATE_estimates,
                                 FI_ATE_estimates=FI_ATE_estimates)
    return(boot_one_time_output)
  }
  
  # Repeat the 'boot_one_time' function B times using parallel computation
  boot_iterations <- function(merged_data, y1_names, y0_names, IV, M_B, para_curr, EM, covariate, B, seed) {
    
    t<-proc.time()
    stopImplicitCluster2 <- function() {
      options <- doParallel:::.options
      if (exists(".revoDoParCluster", where = options) &&
          !is.null(get(".revoDoParCluster", envir = options))) {
        stopCluster(get(".revoDoParCluster", envir = options))
        remove(".revoDoParCluster", envir = options)
      }
    }
    
    registerDoParallel(max(detectCores()-1,1))
    registerDoRNG(seed)
    Sys.time()
    
    ptime <- system.time({
      output1 <- foreach(it = 1:B, .export=ls(envir=parent.frame()), .packages=loadedNamespaces(), .errorhandling = "pass") %dopar% {
        return(boot_one_time(merged_data, y1_names, y0_names, IV, M_B, para_curr, EM, covariate, it=it))
        invisible(gc())
        invisible(gc())
      }
    })[3]
    stopImplicitCluster2()
    iter_time <- proc.time() - t 
    
    # ATT_(b), ATC_(b), ATE_(b) calculated from b-th bootstrap sample
    reg_ATE_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$reg_ATE_estimates))
    FI_ATE_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$FI_ATE_estimates))
    PSM_ATE_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$PSM_ATE_estimates))
    
    # Calculate bootstrap variance estimates
    var_ATE_estimates_reg <- apply(reg_ATE_estimates, MARGIN = 2, FUN = var)
    var_ATE_estimates_FI <- apply(FI_ATE_estimates, MARGIN = 2, FUN = var)
    var_ATE_estimates_PSM <- apply(PSM_ATE_estimates, MARGIN = 2, FUN = var)
    
    boot_iterations_out <- list(PSM_ATE_estimates=PSM_ATE_estimates,
                                reg_ATE_estimates=reg_ATE_estimates,
                                FI_ATE_estimates=FI_ATE_estimates,
                                var_ATE_estimates_reg=var_ATE_estimates_reg,
                                var_ATE_estimates_FI=var_ATE_estimates_FI,
                                var_ATE_estimates_PSM=var_ATE_estimates_PSM,
                                iter_time=iter_time)
    return(boot_iterations_out)
  }
  boot_iterations_out <- boot_iterations(merged_data, y1_names, y0_names, IV, M_B, para_curr, EM, covariate, B, seed)
  return(boot_iterations_out)
}
