# Scale and Manipulate a Matrix

Scales a gene expression matrix, optionally applies logarithmic
transformation, and performs feature manipulation to remove problematic
variables (e.g., genes with zero variance or missing values).

## Usage

``` r
scale_matrix(matrix, log2matrix = TRUE, manipulate = TRUE)
```

## Arguments

- matrix:

  Numeric matrix with genes as rows and samples as columns.

- log2matrix:

  Logical indicating whether to apply log2 transformation using
  \[log2eset()\]. Default is \`TRUE\`.

- manipulate:

  Logical indicating whether to perform feature manipulation to remove
  problematic features. Default is \`TRUE\`.

## Value

A scaled matrix (genes as rows, samples as columns).

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
test_matrix <- matrix(
  rlnorm(100),
  ncol = 5,
  dimnames = list(paste0("Gene", 1:20), paste0("Sample", 1:5))
)
eset2 <- scale_matrix(test_matrix, log2matrix = FALSE, manipulate = TRUE)
#> ✔ Retained 20 of 20 features
#> ℹ Retained 20 features after manipulation
head(eset2)
#>           Sample1    Sample2    Sample3    Sample4    Sample5
#> Gene1 -0.45120149 -0.9503647 -0.6087774  1.5061076  0.5042360
#> Gene2 -0.30999533 -0.2805849 -0.2560728 -0.8807423  1.7273953
#> Gene3  1.77959701 -0.5253800 -0.5654074 -0.3374958 -0.3513138
#> Gene4 -0.40723962 -0.5741448  1.7620081 -0.6084337 -0.1721900
#> Gene5 -0.07783652 -0.5723871  1.7342497 -0.7306081 -0.3534179
#> Gene6  1.72922021 -0.7192169 -0.6553447 -0.1861758 -0.1684828
```
