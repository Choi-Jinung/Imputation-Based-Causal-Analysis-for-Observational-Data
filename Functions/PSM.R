library(MatchIt)
library(marginaleffects)

# Propensity Score Matching (PSM)

# Arguments
## 1. data
## 2. covariates_ps: string vector for names of predictor variables to be included in the propensity score model
## 3. method: a string to specify the matching method
## 4. estimand: average treatment effect to be estimated. "ATC" or "ATT" or "ATE"
## 5. covariates_out: string vector for names of predictor variables to be included in the outcome regression model
## 6. estimand_group: specify the group over which the causal effect will be estimated
## To ensure error-free results, the following rules must be satisfied:
## if estimand=="ATC", estimand_group <- 0
## if estimand=="ATT", estimand_group <- 1
## if estimand=="ATE", estimand_group <- "total"

PSM <- function(data, covariates_ps, method, estimand, covariates_out, estimand_group) {
  formula_ps <- as.formula(paste("t", "~", paste(covariates_ps, collapse = "+")))
  if (method=="nearest") {
    m.out1 <- matchit(formula_ps, data = data, method=method, distance="glm", estimand = estimand)
  } else {
    m.out1 <- matchit(formula_ps, data = data, method=method, distance="glm", estimand = estimand)
  }
  m.data <- match.data(m.out1)
  if (is.null(covariates_out)) {
    formula_out <- as.formula(paste("y", "~", "t"))
  } else {
    formula_out <- as.formula(paste("y", "~", "t", "*", "(", paste(covariates_out, collapse = "+"), ")"))
  }
  out_fit <- lm(formula_out, data=m.data, weights = weights)
  if (estimand_group == "total") {
    est <- avg_comparisons(out_fit,
                           variables = "t",
                           vcov = ~subclass,
                           newdata = m.data,
                           wts = "weights")
  } else {
    est <- avg_comparisons(out_fit,
                           variables = "t",
                           vcov = ~subclass,
                           newdata = subset(m.data, t == estimand_group),
                           wts = "weights")
  }
  output <- list(m.out1=m.out1, m.data=m.data, out_fit=out_fit, est=est)
  return(output)
}