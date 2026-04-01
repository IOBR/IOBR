# Elastic Net Model Fitting

Fits elastic net model with cross-validation to find optimal alpha and
lambda. Searches across a grid of alpha values (0 to 1) and lambda
values to minimize cross-validation error.

## Usage

``` r
Enet(train.x, train.y, lambdamax, nfold = 10)
```

## Arguments

- train.x:

  Training predictors matrix.

- train.y:

  Training binary outcomes (0/1 or factor).

- lambdamax:

  Maximum lambda value for the grid search.

- nfold:

  Number of CV folds. Default is \`10\`.

## Value

List containing:

- chose_alpha:

  Optimal alpha value (0-1)

- chose_lambda:

  Optimal lambda value

## Examples

``` r
# \donttest{
if (requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(123)
  train_data <- matrix(rnorm(50 * 5), ncol = 5)
  train_outcome <- rbinom(50, 1, 0.5)
  result <- Enet(train.x = train_data, train.y = train_outcome, lambdamax = 1, nfold = 5)
}
# }
```
