# Feature Selection for Predictive or Prognostic Models Using LASSO Regression

Applies LASSO (Least Absolute Shrinkage and Selection Operator)
regression to construct predictive or prognostic models. Supports both
binary and survival response variables, utilizing cross-validation for
optimal model selection.

## Usage

``` r
lasso_select(
  x,
  y,
  type = c("binary", "survival"),
  nfold = 10,
  lambda = c("lambda.min", "lambda.1se")
)
```

## Arguments

- x:

  A numeric matrix. Features (e.g., gene symbols or CGI) as row names
  and samples as column names.

- y:

  A response variable vector. Can be binary (0/1) or survival data
  (e.g., survival time and event status).

- type:

  Character. Model type: "binary" for binary response or "survival" for
  survival analysis. Default is "binary".

- nfold:

  Integer. Number of folds for cross-validation. Default is 10.

- lambda:

  Character. Regularization parameter selection: "lambda.min" (minimum
  mean cross-validated error) or "lambda.1se" (one standard error from
  minimum). Default is "lambda.min".

## Value

Character vector of selected feature names with non-zero coefficients in
the optimal LASSO model.

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
gene_expression <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
rownames(gene_expression) <- paste0("Gene", 1:100)
colnames(gene_expression) <- paste0("Sample", 1:20)

# Binary response example
binary_outcome <- sample(c(0, 1), 20, replace = TRUE)
lasso_select(
  x = gene_expression,
  y = binary_outcome,
  type = "binary",
  nfold = 5
)
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> character(0)
```
