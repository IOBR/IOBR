# Calculate Time-Dependent AUC for Survival Models

Evaluates prognostic ability of a survival model by calculating
time-dependent AUC at the 30th and 90th percentiles of survival time.
These thresholds assess short-term and long-term predictive accuracy.

## Usage

``` r
PrognosticAUC(model, newx, s, acture.y)
```

## Arguments

- model:

  A fitted survival model object capable of generating risk scores.

- newx:

  A matrix or data frame of new data for prediction.

- s:

  Lambda value for prediction. Can be numeric or
  \`"lambda.min"\`/\`"lambda.1se"\`.

- acture.y:

  Data frame with \`time\` and \`status\` columns.

## Value

A data frame with AUC values at 30th (\`probs.3\`) and 90th
(\`probs.9\`) percentiles.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
if (requireNamespace("glmnet", quietly = TRUE) &&
  requireNamespace("survival", quietly = TRUE) &&
  requireNamespace("timeROC", quietly = TRUE)) {
  library(survival)
  set.seed(123)
  x <- matrix(rnorm(100 * 5), ncol = 5)
  y <- Surv(rexp(100), rbinom(100, 1, 0.5))
  fit <- glmnet::cv.glmnet(x, y, family = "cox")
  acture_y <- data.frame(time = y[, 1], status = y[, 2])
  auc_results <- PrognosticAUC(fit, newx = x, s = "lambda.min", acture.y = acture_y)
}
# }
```
