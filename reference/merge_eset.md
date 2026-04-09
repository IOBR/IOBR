# Merge Expression Sets by Row Names

Merges two or three expression sets (matrices or data frames) by row
names (gene symbols), removing duplicates. The function ensures common
genes across all input expression sets are retained.

## Usage

``` r
merge_eset(eset1, eset2, eset3 = NULL)
```

## Arguments

- eset1:

  First expression set (matrix or data frame with row names).

- eset2:

  Second expression set (matrix or data frame with row names).

- eset3:

  Optional third expression set. Default is \`NULL\`.

## Value

Merged expression set (data frame) with duplicates removed. Row names
correspond to gene symbols.

## Author

Dongqiang Zeng

## Examples

``` r
# Load example data
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"

# Create mock expression sets with common genes
common_genes <- c("TP53", "BRCA1", "EGFR", "MYC")
eset1 <- matrix(rnorm(12),
  nrow = 4,
  dimnames = list(common_genes, paste0("S", 1:3))
)
eset2 <- matrix(rnorm(16),
  nrow = 4,
  dimnames = list(common_genes, paste0("S", 4:7))
)

# Merge two expression sets
merged_eset <- merge_eset(eset1, eset2)
#> ℹ Common genes between eset1 and eset2: 4
#> ℹ No duplicate gene symbols found.
#> ✔ Final merged expression set: 4 genes x 7 samples
print(dim(merged_eset))
#> [1] 4 7
```
