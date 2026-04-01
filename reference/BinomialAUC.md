# Calculate AUC for Binomial Model

Computes Area Under the ROC Curve (AUC) for model predictions using the
ROCR package. Handles binary classification models from glmnet.

## Usage

``` r
BinomialAUC(model, newx, s, acture.y)
```

## Arguments

- model:

  Fitted glmnet model object.

- newx:

  New data matrix for prediction.

- s:

  Lambda value for prediction (e.g., "lambda.min" or numeric).

- acture.y:

  Actual binary outcomes (numeric 0/1 or factor).

## Value

Numeric AUC value between 0 and 1.

## Examples

``` r
# \donttest{
if (requireNamespace("glmnet", quietly = TRUE) && requireNamespace("ROCR", quietly = TRUE)) {
  set.seed(123)
  train_data <- matrix(rnorm(100 * 5), ncol = 5)
  train_outcome <- rbinom(100, 1, 0.5)
  test_data <- matrix(rnorm(50 * 5), ncol = 5)
  test_outcome <- rbinom(50, 1, 0.5)
  fitted_model <- glmnet::cv.glmnet(train_data, train_outcome, family = "binomial", nfolds = 5)
  auc_value <- BinomialAUC(fitted_model, test_data, fitted_model$lambda.min, test_outcome)
  print(auc_value)
}
#> [1] 0.4879227
# }
```
