getwd()
setwd(getwd())
rm(list=ls())

# Import functions
# source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/fractional imputation_single.R")
# source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/variance_estimation_single.R")
source("/Users/a82103/Desktop/Research/Causal Inference/Code_final/Functions/fractional imputation_single.R")
source("/Users/a82103/Desktop/Research/Causal Inference/Code_final/Functions/variance_estimation_single.R")

# Import the control group and the treatment group data 
data_t0 <- read.csv("data_t0.csv")[, c("y0", "x1", "x2", "x3", "x4")] %>% rename(y=y0)
data_t1 <- read.csv("data_t1.csv")[, c("y1", "x1", "x2", "x3", "x4")] %>% rename(y=y1)

# Function to generate real data analysis outcomes
## Arguments
### 1. data_t0: control group data
### 2. data_t1: treatment group data
### 3. IV: instrumental variables
### 4. M_B: imputation size
real.data.analysis <- function(data_t0, data_t1, IV, M_B) {
  
  # Impute the control group data, data_t0
  # dataB: data that will be fractionally imputed
  dataA <- data_t1; dataB <- data_t0
  dataA <- dataA %>% rename(y1=y)
  dataB <- dataB %>% rename(y2=y)
  p <- ncol(dataA)
  beta_init <- rep(0, p-length(IV)+1)
  var_init <- 1
  # Apply fractional imputation to the control group data
  imputed.data_t0_list <- FI_single(dataA, dataB, IV, M_B, c(beta_init, var_init), EM=T)
  # estimated parameter of f(y0 | x, y1)
  est_para_t0 <- imputed.data_t0_list$est_para 
  # 95% confidence intervals of the above estimated parameters
  CI_y0_xy1 <- variance_estimation_single(imputed.data_t0_list)$CI_95 
  # estimated missing potential outcome under treatment
  y1.hat <- imputed.data_t0_list$imputed_data %>% group_by(id_B) %>% summarise(y1.hat = sum(w*y1)) %>% .[[2]] 
  # Average Treatment Effect under Control (ATC) 
  ATE_t0 <- mean(y1.hat) - mean(data_t0$y)
  
  # Impute the treatment group data, data_t1
  dataA <- data_t0; dataB <- data_t1
  dataA <- dataA %>% rename(y1=y)
  dataB <- dataB %>% rename(y2=y)
  p <- ncol(dataA)
  beta_init <- rep(0, p-length(IV)+1)
  var_init <- 1
  # Apply fractional imputation to the treatment group data
  imputed.data_t1_list <- FI_single(dataA, dataB, IV, M_B, c(beta_init, var_init), EM=T)
  # estimated parameter of f(y1 | x, y0)
  est_para_t1 <- imputed.data_t1_list$est_para
  # 95% confidence intervals of the above estimated parameters
  CI_y1_xy0 <- variance_estimation_single(imputed.data_t1_list)$CI_95
  # estimated missing potential outcome under control
  y0.hat <- imputed.data_t1_list$imputed_data %>% group_by(id_B) %>% summarise(y0.hat = sum(w*y1)) %>% .[[2]]
  # Average Treatment Effect under Treatment (ATT) 
  ATE_t1 <- mean(data_t1$y) - mean(y0.hat)
  
  # Average Treatment Effect (ATE)
  ATE_total <- mean(c(data_t1$y, y1.hat)-c(y0.hat, data_t0$y))
    
  result <- list(ATE_t0=ATE_t0, ATE_t1=ATE_t1, ATE_total=ATE_total, 
                 est_para_t0=est_para_t0, est_para_t1=est_para_t1,
                 CI_y0_xy1=CI_y0_xy1, CI_y1_xy0=CI_y1_xy0,
                 imputed.data_t0_list=imputed.data_t0_list, imputed.data_t1_list=imputed.data_t1_list)
  return(result)
}

# Get the real data analysis outcome for each sensitivity analysis scenario
set.seed(20240805)
est_result1 <- real.data.analysis(data_t0 = data_t0, data_t1 = data_t1,
                                  IV = c("x1", "x2", "x3", "x4"), 
                                  M_B = 2000)
est_result2 <- real.data.analysis(data_t0 = data_t0, data_t1 = data_t1,
                                  IV = c("x2", "x3", "x4"), 
                                  M_B = 2000)
est_result3 <- real.data.analysis(data_t0 = data_t0, data_t1 = data_t1,
                                  IV = c("x3", "x4"), 
                                  M_B = 2000)

# Summarize ATT, ATC, ATE estimates from FICI by scenario
result_ATE <- rbind(c(est_result1$ATE_t1, est_result1$ATE_t0, est_result1$ATE_total),
                    c(est_result2$ATE_t1, est_result2$ATE_t0, est_result2$ATE_total),
                    c(est_result3$ATE_t1, est_result3$ATE_t0, est_result3$ATE_total))
colnames(result_ATE) <- c("ATT", "ATC", "ATE")
rownames(result_ATE) <- c("scenario1", "scenario2", "scenario3")
round(result_ATE, digits = 2)

# 95% confidence intervals for the parameters of f(y0 | x, y1)
round(est_result1$CI_y0_xy1, digits = 2)
round(est_result2$CI_y0_xy1, digits = 2)
round(est_result3$CI_y0_xy1, digits = 2)

# 95% confidence intervals for the parameters of f(y1 | x, y0)
round(est_result1$CI_y1_xy0, digits = 2)
round(est_result2$CI_y1_xy0, digits = 2)
round(est_result3$CI_y1_xy0, digits = 2)
