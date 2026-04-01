# Differential Expression Analysis Using Limma

Performs differential expression analysis using the limma package on a
given gene expression dataset. Constructs a design matrix from phenotype
data, fits a linear model, applies contrasts, and computes statistics
for differential expression.

## Usage

``` r
limma.dif(exprdata, pdata, contrastfml)
```

## Arguments

- exprdata:

  A matrix with rownames as features like gene symbols or cgi, and
  colnames as samples.

- pdata:

  A two-column dataframe where the first column matches the colnames of
  exprdata and the second column contains the grouping variable.

- contrastfml:

  A character vector for contrasts to be tested (see ?makeContrasts for
  more details).

## Value

Returns a dataframe from limma::topTable, which includes genes as rows
and columns like genelist, logFC, AveExpr, etc.

## Examples

``` r
# Toy example with 100 genes and 6 samples
set.seed(123)
exprdata <- matrix(
  rnorm(100 * 6),
  nrow = 100,
  ncol = 6,
  dimnames = list(
    paste0("gene", 1:100),
    paste0("sample", 1:6)
  )
)

# Phenotype data: 3 vs 3
pdata <- data.frame(
  sample = colnames(exprdata),
  group = rep(c("group1", "group2"), each = 3),
  stringsAsFactors = FALSE
)

# Differential expression: group1 vs group2
res <- limma.dif(
  exprdata = exprdata,
  pdata = pdata,
  contrastfml = "group1 - group2"
)
head(res)
#>            logFC     AveExpr         t    P.Value adj.P.Val         B
#> gene16  2.146063  0.18155652  2.637222 0.01296770 0.6279722 -4.554407
#> gene91  1.979646 -0.19771191  2.533719 0.01657904 0.6279722 -4.557543
#> gene97  1.791447  0.65590789  2.412055 0.02199743 0.6279722 -4.561178
#> gene71 -1.738285  0.91534298 -2.330224 0.02650410 0.6279722 -4.563586
#> gene88  1.551421 -0.05545564  2.105013 0.04352443 0.6279722 -4.570024
#> gene56  1.644500 -0.43974674  2.038109 0.05018179 0.6279722 -4.571876
```
