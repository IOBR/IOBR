# Convert Mouse Gene Symbols to Human Gene Symbols

Converts mouse gene symbols to human gene symbols in an expression
dataset. Supports using either an online resource (Ensembl) or a local
dataset for conversion.

## Usage

``` r
mouse2human_eset(
  eset,
  source = c("local", "ensembl"),
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
if (interactive()) {
  set.seed(123)
  data <- matrix(runif(50 * 3), nrow = 50, ncol = 3)
  rownames(data) <- c("Tpt1", "Hmgb1", "Gapdh", paste0("Gene", 4:50))
  colnames(data) <- paste0("Sample", 1:3)
  human_data <- mouse2human_eset(data, source = "local", is_matrix = TRUE)
  if (!is.null(human_data)) head(human_data)
}
```
