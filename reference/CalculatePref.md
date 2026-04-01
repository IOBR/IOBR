# Calculate Performance Metrics

Computes True Positive Rate (TPR) and False Positive Rate (FPR) for ROC
analysis using the ROCR package. Used internally for ROC curve
generation.

## Usage

``` r
CalculatePref(model, newx, s, acture.y)
```

## Arguments

- model:

  Fitted glmnet model.

- newx:

  New data matrix for prediction.

- s:

  Lambda value for prediction.

- acture.y:

  Actual binary outcomes.

## Value

ROCR performance object containing TPR and FPR values.

## Examples

``` r
# \donttest{
if (requireNamespace("glmnet", quietly = TRUE) && requireNamespace("ROCR", quietly = TRUE)) {
  fitted_model <- glmnet::cv.glmnet(matrix(rnorm(100), ncol = 2), rbinom(50, 1, 0.5), nfolds = 3)
  perf <- CalculatePref(fitted_model, matrix(rnorm(20), ncol = 2), "lambda.min", rbinom(10, 1, 0.5))
}
# }
```
