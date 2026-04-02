# Compute Prognostic Results for Survival Models

Computes and compiles prognostic results from a survival model fitted
with \`glmnet\`. Extracts model coefficients at optimal lambda values
(\`lambda.min\` and \`lambda.1se\`) and calculates time-dependent AUC
metrics for both training and testing datasets.

## Usage

``` r
PrognosticResult(model, train.x, train.y, test.x, test.y)
```

## Arguments

- model:

  A fitted survival model object (e.g., from \`glmnet::cv.glmnet\`).

- train.x:

  Matrix or data frame of training predictors.

- train.y:

  Training dataset survival outcomes (time and status).

- test.x:

  Matrix or data frame of testing predictors.

- test.y:

  Testing dataset survival outcomes (time and status).

## Value

A list containing:

- model:

  The fitted model object

- coefs:

  Data frame of coefficients at \`lambda.min\` and \`lambda.1se\`

- AUC:

  Data frame with AUC values for train/test at both lambda values

## Author

Dongqiang Zeng

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE) &&
  requireNamespace("survival", quietly = TRUE) &&
  requireNamespace("timeROC", quietly = TRUE)) {
  library(survival)
  set.seed(123)
  train_x <- matrix(rnorm(100 * 10), ncol = 10)
  train_y <- data.frame(time = rexp(100), status = rbinom(100, 1, 0.5))
  test_x <- matrix(rnorm(50 * 10), ncol = 10)
  test_y <- data.frame(time = rexp(50), status = rbinom(50, 1, 0.5))
  fit <- glmnet::cv.glmnet(train_x, Surv(train_y$time, train_y$status), family = "cox")
  results <- PrognosticResult(
    model = fit, train.x = train_x, train.y = train_y,
    test.x = test_x, test.y = test_y
  )
}
```
