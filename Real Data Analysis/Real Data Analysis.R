getwd()
setwd(getwd())
rm(list=ls())
library(devtools)
install.packages("/Users/a82103/Desktop/MatchIt", repos = NULL, type = "source")
library(MatchIt)

# Import functions
source("/Users/choi-jinung/Desktop/Research/Causal Inference/Code_final/Functions/bootstrap.R")
# source("/Users/a82103/Desktop/Research/Causal Inference/Code_final/Functions/bootstrap.R")

# Data description
data_t0 <- read.csv("data_t0.csv")[, c("y0", "x1", "x2", "x3", "x4")]
data_t1 <- read.csv("data_t1.csv")[, c("y1", "x1", "x2", "x3", "x4")]
n1 <- nrow(data_t1); n0 <- nrow(data_t0)
summary(data_t1)
summary(data_t0)
back_data_t1 <- exp(data_t1)-1
back_data_t0 <- exp(data_t0)-1
summary(back_data_t1)
summary(back_data_t0)

# Summarize the matching outcome
data_t0_matching <- data_t0 %>% mutate(t=0) %>% rename(y=y0, DURA=x1, IMG=x2, VDO=x3, LENGTH=x4)
data_t1_matching <- data_t1 %>% mutate(t=1) %>% rename(y=y1, DURA=x1, IMG=x2, VDO=x3, LENGTH=x4)
data_matching <- rbind(data_t0_matching, data_t1_matching)
m.out <- matchit(t~DURA+IMG+VDO+LENGTH, data = data_matching, method="full", distance="glm", estimand = "ATE")
summary(m.out)
plot(m.out, type = "density", interactive = F, which.xs = ~DURA + IMG + VDO)
plot(m.out, type = "density", interactive = F, which.xs = ~LENGTH)

# Diagnosis of the fitted proposal distribution
lm_treatment <- lm(y1~x1+x2+x3+x4, data=data_t1)
lm_control <- lm(y0~x1+x2+x3+x4, data=data_t0)
summary(lm_treatment)
summary(lm_control)
par(mfrow=c(1,2))
plot(lm_treatment)
plot(lm_control)

# Variance estimation
var_est_out <- var_est(control = data_t0, treatment = data_t1,
                       IV = c("x2", "x3", "x4"),
                       M_B=100, para_curr = c(0, 0, 0, 1), EM=T, 
                       covariate = c("x1", "x2", "x3", "x4"),
                       B=1000, seed=10)

# Point estimation
# CIA
reg_ATT_list <- Regression_CI(data_t0, data_t1, c("x1", "x2", "x3", "x4"))
reg_ATC_list <- Regression_CI(data_t1, data_t0, c("x1", "x2", "x3", "x4"))
reg_ATT <- reg_ATT_list$observed_bar - reg_ATT_list$missing_bar
reg_ATC <- reg_ATC_list$missing_bar - reg_ATC_list$observed_bar
reg_ATE <- (n1*reg_ATT + n0*reg_ATC)/(n1+n0)
reg_ATE_estimates <- c(reg_ATT, reg_ATC, reg_ATE)

# PSM
control_psm <- data.frame(data_t0) %>% rename(y=y0) %>% mutate(t=0)
treatment_psm <- data.frame(data_t1) %>% rename(y=y1) %>% mutate(t=1)
data_psm <- rbind(control_psm, treatment_psm)
PSM_ATT <- PSM(data=data_psm,
               method="full", estimand = "ATT", estimand_group = 1,
               covariates_ps = c("x1","x2","x3","x4"),
               covariates_out = c("x1","x2","x3","x4"))
PSM_ATC <- PSM(data=data_psm,
               method="full", estimand = "ATC", estimand_group = 0,
               covariates_ps = c("x1","x2","x3","x4"),
               covariates_out = c("x1","x2","x3","x4"))
PSM_ATE <- PSM(data=data_psm,
               method="full", estimand = "ATE", estimand_group = "total",
               covariates_ps = c("x1","x2","x3","x4"),
               covariates_out = c("x1","x2","x3","x4"))
PSM_ATE_estimates <- c(PSM_ATT$est$estimate, PSM_ATC$est$estimate, PSM_ATE$est$estimate)

# FI
set.seed(1)
impute_control <- FI_single(dataA=data_t1, dataB=data_t0, IV=c("x2", "x3", "x4"), M_B=2000, para_curr=c(0,0,0,1), EM=T)
impute_treatment <- FI_single(dataA=data_t0, dataB=data_t1, IV=c("x2", "x3", "x4"), M_B=2000, para_curr=c(0,0,0,1), EM=T)
FI_ATC <- impute_control$missing_bar - impute_control$observed_bar
FI_ATT <- impute_treatment$observed_bar - impute_treatment$missing_bar
FI_ATE <- (n0*FI_ATC + n1*FI_ATT)/(n0+n1)
FI_ATE_estimates <- c(FI_ATT, FI_ATC, FI_ATE)

# Calculate p-values
p_values_reg <- (1-pnorm(q=reg_ATE_estimates, mean=0, sd=sqrt(var_est_out$var_ATE_estimates_reg)))*2
p_values_PSM <- (1-pnorm(q=PSM_ATE_estimates, mean=0, sd=sqrt(var_est_out$var_ATE_estimates_PSM)))*2
p_values_FI <- (1-pnorm(q=FI_ATE_estimates, mean=0, sd=sqrt(var_est_out$var_ATE_estimates_FI)))*2

# Summarize the result
ATE_summary <- round(rbind(reg_ATE_estimates, PSM_ATE_estimates, FI_ATE_estimates), digits = 2)
rownames(ATE_summary) <- c("CIA", "PSM", "FI")
colnames(ATE_summary) <- c("ATT", "ATC", "ATE")

p_val_summary <- round(rbind(p_values_reg, p_values_PSM, p_values_FI), digits = 2)
rownames(p_val_summary) <- c("CIA", "PSM", "FI")
colnames(p_val_summary) <- c("ATT", "ATC", "ATE")
