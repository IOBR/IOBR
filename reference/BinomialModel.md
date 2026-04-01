# Binomial Model Construction

Constructs and evaluates binomial logistic regression models using Lasso
and Ridge regularization. Processes input data, scales features if
specified, splits data into training/testing sets, and fits both Lasso
and Ridge models. Optionally generates AUC plots for model evaluation.

## Usage

``` r
BinomialModel(
  x,
  y,
  seed = 123456,
  scale = TRUE,
  train_ratio = 0.7,
  nfold = 10,
  plot = TRUE,
  palette = "jama",
  cols = NULL
)
```

## Arguments

- x:

  A data frame containing sample ID and features. First column must be
  sample ID.

- y:

  A data frame where first column is sample ID and second column is
  outcome (numeric or factor).

- seed:

  Integer for random seed. Default is \`123456\`.

- scale:

  Logical indicating whether to scale features. Default is \`TRUE\`.

- train_ratio:

  Numeric between 0 and 1 for training proportion. Default is \`0.7\`.

- nfold:

  Integer for cross-validation folds. Default is \`10\`.

- plot:

  Logical indicating whether to generate AUC plots. Default is \`TRUE\`.

- palette:

  Character string for color palette. Default is \`"jama"\`.

- cols:

  Optional color vector for ROC curves. Default is \`NULL\`.

## Value

List containing:

- lasso_result:

  Lasso model results

- ridge_result:

  Ridge model results

- train.x:

  Training data with IDs

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
imvigor210_sig <- load_data("imvigor210_sig")
imvigor210_pdata <- load_data("imvigor210_pdata")
pdata_group <- imvigor210_pdata[imvigor210_pdata$BOR_binary != "NA", c("ID", "BOR_binary")]
pdata_group$BOR_binary <- factor(ifelse(pdata_group$BOR_binary == "R", 1, 0))
result <- BinomialModel(
  x = imvigor210_sig, y = pdata_group,
  seed = 123456, scale = TRUE, train_ratio = 0.7, nfold = 10, plot = FALSE
)
#> ℹ Processing data
#> ℹ Splitting data into training and test sets
#> ℹ Running LASSO
#> ℹ Running RIDGE REGRESSION
#> ✔ Model fitting complete
# }
```
