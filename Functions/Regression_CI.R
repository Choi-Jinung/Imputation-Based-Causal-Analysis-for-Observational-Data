# Regression Imputation Estimator under Conditional Independence Assumption (CIA)

# Arguments
## 1. dataA: data frame with X and "y1"
## 2. dataB: data frame with X and "y2". To predict y1 values of this data set is the goal.
## 3. covariate: string vector for names of predictor variables to be included in the regression model.

Regression_CI <- function(dataA, dataB, covariate) {
  y1_names <- setdiff(colnames(dataA), colnames(dataB))
  y2_names <- setdiff(colnames(dataB), colnames(dataA))
  formula_reg <- as.formula( paste(y1_names, "~", paste(covariate, collapse = "+")) )
  prop_fit <- lm(formula_reg, data=dataA)
  missing_hat <- predict(prop_fit, newdata=dataB)
  missing_bar <- mean(missing_hat)
  observed_bar <- mean(dataB[, y2_names])
  out <- list(prop_fit=prop_fit, missing_hat=missing_hat, missing_bar=missing_bar, observed_bar=observed_bar)
  return(out)
}