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
# \donttest{
# Load and annotate expression data
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
anno_rnaseq <- load_data("anno_rnaseq")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "anno_rnaseq"
eset_stad <- anno_eset(eset = eset_stad, annotation = anno_rnaseq)
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2098 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50139
eset_stad <- tibble::rownames_to_column(as.data.frame(eset_stad), var = "id")

# Create duplicate gene names for demonstration
eset_stad[2:3, "id"] <- "MT-CO1"

# Check duplicates before
sum(duplicated(eset_stad$id))
#> [1] 2

# Remove duplicates using mean expression as ranking criterion
eset_stad <- remove_duplicate_genes(
  eset = eset_stad,
  column_of_symbol = "id",
  method = "mean"
)
#> ℹ Found 2 duplicate symbols. Using "mean" for ranking.
#> ✔ Reduced to 50137 unique genes

# Check duplicates after
sum(duplicated(rownames(eset_stad)))
#> [1] 0
# }
```
