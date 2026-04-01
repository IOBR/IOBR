# Regression Result Computation

Computes regression results with coefficients at lambda.min and
lambda.1se, and evaluates AUC for binomial outcomes. Returns a
comprehensive summary of model performance on both training and testing
datasets.

## Usage

``` r
RegressionResult(train.x, train.y, test.x, test.y, model)
```

## Arguments

- train.x:

  Training predictors matrix.

- train.y:

  Training outcomes (binary factor).

- test.x:

  Testing predictors matrix.

- test.y:

  Testing outcomes (binary factor).

- model:

  Fitted cv.glmnet model object.

## Value

List containing:

- model:

  The fitted cv.glmnet model

- coefs:

  Data frame with feature names and coefficients at lambda.min and
  lambda.1se

- AUC:

  Matrix of AUC values for train/test sets at both lambda values

## Examples

``` r
# \donttest{
if (requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(123)
  train_data <- matrix(rnorm(100 * 10), ncol = 10)
  train_outcome <- rbinom(100, 1, 0.5)
  test_data <- matrix(rnorm(50 * 10), ncol = 10)
  test_outcome <- rbinom(50, 1, 0.5)
  fitted_model <- glmnet::cv.glmnet(train_data, train_outcome, family = "binomial", nfolds = 5)
  results <- RegressionResult(
    train.x = train_data, train.y = train_outcome,
    test.x = test_data, test.y = test_outcome, model = fitted_model
  )
}
# }
```
