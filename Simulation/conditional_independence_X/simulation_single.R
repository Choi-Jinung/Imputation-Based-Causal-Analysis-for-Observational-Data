setwd(getwd())
rm(list=ls())
source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/fractional imputation_single.R")
source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/Regression_CI.R")
source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/PSM.R")

library('doParallel')
library('doRNG')
library(MASS)

simulation <- function(nothing) {
  # Generate data
  n <- 800
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y0 <- 1 + x1 + x2 + rnorm(n)
  XX_out <- cbind(1, x1, y0)
  beta0 <- 1; beta1 <- 1; beta2 <- 1; beta <- c(beta0, beta1, beta2); sig <- 1
  y1 <- XX_out%*%beta + rnorm(n, sd=sig)
  full_sample <- data.frame(x1 = x1, x2 = x2, y0 = y0, y1 = y1)
  
  # Strongly ignorable treatment assignment
  alpha0 <- 0; alpha1 <- 1; alpha2 <- 1; alpha <- c(alpha0, alpha1, alpha2)
  XX_res <- cbind(1, x1, x2)
  inclusion_prob <- exp(XX_res%*%alpha)/(1+exp(XX_res%*%alpha))
  r <- rbinom(n=n, size=1, prob = inclusion_prob)
  
  control <- full_sample[r==0, ]; control$y1 <- NULL # control group
  treatment <- full_sample[r==1, ]; treatment$y0 <- NULL # treatment group
  
  # Regression Imputation (CIA)
  reg_out <- Regression_CI(dataA=control, dataB=treatment, covariate=c("x1", "x2"))
  MSE_reg <- mean((full_sample[r==1, ]$y0 - reg_out$missing_hat)^2)
  MAE_reg <- mean(abs(full_sample[r==1, ]$y0 - reg_out$missing_hat))
  treatment_reg <- data.frame(treatment, y0=reg_out$missing_hat)
  reg_ATE <- reg_out$observed_bar - reg_out$missing_bar
  
  # Full Matching PSM Imputation (PSM)
  treatment_PSM <- treatment %>% rename(y=y1) %>% mutate(t=1)
  control_PSM <- control %>% rename(y=y0) %>% mutate(t=0)
  data_PSM <- rbind(treatment_PSM, control_PSM)
  psm_out_full <- PSM(data=data_PSM, covariates_ps=c("x1", "x2"), method="full", estimand="ATT", covariates_out=c("x1", "x2"), estimand_group=1)
  treatment_full <- psm_out_full$m.data %>% filter(t==1)
  treatment0_full <- treatment_full %>% mutate(t=0)
  y1_hat <- predict(psm_out_full$out_fit, newdata = treatment_full)
  y0_hat <- predict(psm_out_full$out_fit, newdata = treatment0_full)
  ITE_PSM <- y1_hat - y0_hat
  ITE_true <- full_sample[r==1, ]$y1 - full_sample[r==1, ]$y0 
  MSE_psm_full <- mean((ITE_PSM - ITE_true)^2)
  MAE_psm_full <- mean(abs(ITE_PSM - ITE_true))
  psm_full_ATE <- psm_out_full$est$estimate
  
  # Fractional Imputation (FI)
  IV <- c("x2")
  M_B <- 100
  beta1_curr <- c(0,0,0); var1_curr <- 1
  para_curr <- c(beta1_curr, var1_curr)
  full_fit <- lm(y1~y0+x1, data=full_sample[r==1, ]); true_var <- sum(full_fit$residuals^2)/full_fit$df.residual
  fi_out <- FI_single(dataA=control, dataB=treatment, IV=IV, M_B=M_B, para_curr=para_curr, EM=T)
  MSE_FI <- mean((full_sample[r==1, ]$y0 - fi_out$missing_hat)^2)
  MAE_FI <- mean(abs(full_sample[r==1, ]$y0 - fi_out$missing_hat))
  full_para <- c(full_fit$coefficients, true_var)
  est_para <- fi_out$est_para
  fi_ATE <- fi_out$observed_bar - fi_out$missing_bar
  true_ATE <- mean(full_sample[r==1, ]$y1 - full_sample[r==1, ]$y0)
  
  MSE <- c(MSE_reg, MSE_psm_full, MSE_FI)
  MAE <- c(MAE_reg, MAE_psm_full, MAE_FI)
  ATE_estimates <- c(true_ATE, reg_ATE, psm_full_ATE, fi_ATE)
  
  sim.result <- list(est_para=est_para, MSE=MSE, MAE=MAE, ATE_estimates=ATE_estimates)
  return(sim.result)
}

# MCMC simulation using parallel computing
B <- 5000

# parallel computation
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
registerDoRNG(2)
Sys.time()

ptime <- system.time({
  output1 <- foreach(it = 1:B, .export=ls(envir=parent.frame()), .packages=loadedNamespaces(), .errorhandling = "pass") %dopar% {
    return(simulation(it))
    invisible(gc())
    invisible(gc())
  }
})[3]
stopImplicitCluster2()
proc.time() - t 

# Summarize the simulation results
est_para_result <- do.call(rbind, lapply(output1, FUN = function(x) x$est_para))
MSE_result <- do.call(rbind, lapply(output1, FUN = function(x) x$MSE))
MAE_result <- do.call(rbind, lapply(output1, FUN = function(x) x$MAE))
ATE_estimates_result <- do.call(rbind, lapply(output1, FUN = function(x) x$ATE_estimates))
colnames(MSE_result) <- c("Reg", "PSM", "FICI")
colnames(MAE_result) <- colnames(MSE_result)
colnames(ATE_estimates_result) <- c("True", colnames(MSE_result))
apply(est_para_result, MARGIN = 2, FUN = mean)

round(apply(ATE_estimates_result, MARGIN = 2, FUN = mean), digits = 2)
round(apply(ATE_estimates_result, MARGIN = 2, FUN = var), digits = 3)
round(apply(MSE_result, MARGIN = 2, FUN = mean), digits = 2)
round(apply(MAE_result, MARGIN = 2, FUN = mean), digits = 2)
par(mfrow=c(1,2))
boxplot(MSE_result, main="MSE")
boxplot(MAE_result, main="MAE")
