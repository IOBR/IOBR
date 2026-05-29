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
# Create example mouse expression data
set.seed(123)
data <- matrix(runif(100 * 5), nrow = 100, ncol = 5)
rownames(data) <- c("Tpt1", "Hmgb1", "Gapdh", paste0("Gene", 4:100))
colnames(data) <- paste0("Sample", 1:5)

# Convert using local database (may return NULL if no internet for internal data)
human_data <- mouse2human_eset(data, source = "local", is_matrix = TRUE)
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "mus_human_gene_symbol"
#> ℹ Row number of original eset: 100
#> ✔ 3% of probes in expression set were annotated
#> ℹ Row number after filtering duplicated gene symbol: 3
if (!is.null(human_data)) head(human_data)
#>         Sample1   Sample2   Sample3     Sample4   Sample5
#> GAPDH 0.4089769 0.4886130 0.6013657 0.779065883 0.9053096
#> HMGB1 0.7883051 0.3328235 0.9623589 0.009429905 0.1370675
#> TPT1  0.2875775 0.5999890 0.2387260 0.784575267 0.9860543
```
