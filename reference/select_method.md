# Select a Signature Scoring Method Subset

Filters an integrated signature score matrix to retain results from a
specified method (PCA, ssGSEA, or zscore) and strips method suffixes
from column names.

## Usage

``` r
select_method(data, method = c("ssGSEA", "PCA", "zscore"))
```

## Arguments

- data:

  Data frame or matrix. Integrated signature score matrix.

- method:

  Character. One of "PCA", "ssGSEA", or "zscore" (case-insensitive).
  Default is "ssGSEA".

## Value

Matrix or data frame containing only the selected method's scores.

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
anno_grch38 <- load_data("anno_grch38")
#> ℹ Loading cached data: "anno_grch38"
hallmark <- load_data("hallmark")
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "hallmark"
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
eset <- eset[1:5000, 1:10]
res <- calculate_sig_score(
  eset = eset,
  signature = hallmark[1:4],
  method = "integration"
)
#> ℹ Calculating signature scores using PCA, z-score, and ssGSEA methods
#> ✔ Applied log2 transformation
#> ℹ Step 1/3: PCA method
#> ℹ Step 2/3: z-score method
#> ℹ Step 3/3: ssGSEA method
#> ℹ GSVA version 2.5.41
#> ℹ Searching for rows with constant values
#> ℹ Calculating ssGSEA scores for 4 gene sets
#> ℹ Normalizing ssGSEA scores
#> ✔ Calculations finished
select_method(res, method = "PCA")
#> # A tibble: 10 × 6
#>    ID           Index HALLMARK_ADIPOGENESIS HALLMARK_ALLOGRAFT_REJECTION
#>    <chr>        <int>                 <dbl>                        <dbl>
#>  1 TCGA-BR-6455     1                 7.41                        6.81  
#>  2 TCGA-BR-7196     2                -1.53                        4.46  
#>  3 TCGA-BR-8371     3                -3.68                       -6.29  
#>  4 TCGA-BR-8380     4                 2.73                       -3.10  
#>  5 TCGA-BR-8592     5                 5.04                       -0.0296
#>  6 TCGA-BR-8686     6               -17.5                        -5.63  
#>  7 TCGA-BR-A4IV     7                 3.73                       -2.48  
#>  8 TCGA-BR-A4J4     8                 4.54                        4.35  
#>  9 TCGA-BR-A4J9     9                -0.214                      -5.36  
#> 10 TCGA-FP-7916    10                -0.501                       7.27  
#> # ℹ 2 more variables: HALLMARK_ANDROGEN_RESPONSE <dbl>,
#> #   HALLMARK_ANGIOGENESIS <dbl>
# }
```
