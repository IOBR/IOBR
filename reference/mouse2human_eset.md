# Convert Mouse Gene Symbols to Human Gene Symbols

Converts mouse gene symbols to human gene symbols in an expression
dataset. Supports using either an online resource (Ensembl) or a local
dataset for conversion.

## Usage

``` r
mouse2human_eset(
  eset,
  source = c("ensembl", "local"),
  is_matrix = TRUE,
  column_of_symbol = NULL,
  verbose = FALSE
)
```

## Arguments

- eset:

  Matrix or data frame. Expression matrix with genes in rows.

- source:

  Character. Data source for conversion: "ensembl" (online) or "local".
  Default is "ensembl". If Ensembl fails, use "local" which uses the
  internal \`mus_human_gene_symbol\` dataset.

- is_matrix:

  Logical. Whether \`eset\` is a matrix with gene symbols as row names.
  Default is \`TRUE\`. If \`FALSE\`, \`column_of_symbol\` must be
  specified.

- column_of_symbol:

  Character or \`NULL\`. Column name containing gene symbols if \`eset\`
  is not a matrix. Default is \`NULL\`.

- verbose:

  Logical. If \`TRUE\`, prints available Ensembl datasets. Default is
  \`FALSE\`.

## Value

Expression set with human gene symbols.

## Author

Dongqiang Zeng

## Examples

``` r
# Create example mouse expression data
anno_gc_vm32 <- load_data("anno_gc_vm32")
num_rows <- 200
num_cols <- 10
sample_names <- paste0("Sample", 1:num_cols)
data <- matrix(runif(num_rows * num_cols), nrow = num_rows, ncol = num_cols)
rownames(data) <- anno_gc_vm32$symbol[1:200]
colnames(data) <- sample_names

# Convert using local database
human_data <- mouse2human_eset(data, source = "local", is_matrix = TRUE)
#> ℹ Row number of original eset: 200
#> ✔ 94% of probes in expression set were annotated
#> ℹ Found 1 duplicate symbol, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 217
```
