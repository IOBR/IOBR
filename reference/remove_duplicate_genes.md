# Remove Duplicate Gene Symbols in Gene Expression Data

This function addresses duplicate gene symbols in a gene expression
dataset by selecting the highest-expressing instance among duplicates.
Users can choose between mean, standard deviation, or sum as the ranking
criterion for selection. This is useful for preparing data where
duplicates can lead to issues in downstream analyses.

## Usage

``` r
remove_duplicate_genes(eset, column_of_symbol, method = c("mean", "sd", "sum"))
```

## Arguments

- eset:

  A data frame or matrix representing gene expression data, with gene
  symbols as one of the columns.

- column_of_symbol:

  The name of the column containing gene symbols in \`eset\`.

- method:

  The ranking method to use for selecting among duplicate gene symbols:
  \`"mean"\` for mean expression, \`"sd"\` for standard deviation, or
  \`"sum"\` for sum of expression values. Default is \`"mean"\`.

## Value

A modified version of \`eset\` where duplicate gene symbols have been
reduced to a single entry (the highest-ranking one). The gene symbols
are set as row names in the returned data frame.

## Note

Important: This function performs selection, not aggregation. For
duplicate genes, it retains only the highest-ranking instance (based on
the specified method) and discards others.

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
test_eset <- data.frame(
  symbol = c("GeneA", "GeneA", "GeneB", "GeneC"),
  S1 = c(10, 5, 20, 15),
  S2 = c(12, 7, 22, 17)
)
# Remove duplicates using mean expression
test_eset_unique <- remove_duplicate_genes(
  eset = test_eset,
  column_of_symbol = "symbol",
  method = "mean"
)
#> ℹ Found 1 duplicate symbol. Using "mean" for ranking.
#> ✔ Reduced to 3 unique genes
print(test_eset_unique)
#>       S1 S2
#> GeneB 20 22
#> GeneC 15 17
#> GeneA 10 12
```
