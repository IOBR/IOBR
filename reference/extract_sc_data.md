# Extract Data Frame from Seurat Object

Extracts and combines a data frame with cells as rows and features as
columns from Seurat assay data. Supports multiple assays and optional
metadata integration.

## Usage

``` r
extract_sc_data(
  sce,
  vars = NULL,
  assay,
  slot = "scale.data",
  combine_meta_data = TRUE
)
```

## Arguments

- sce:

  Seurat object.

- vars:

  Character vector of feature names to extract. If \`NULL\`, all
  features are extracted.

- assay:

  Character vector specifying assay(s) to pull data from.

- slot:

  Character string specifying the assay data slot. Default is
  \`"scale.data"\`.

- combine_meta_data:

  Logical indicating whether to combine metadata with the extracted data
  frame. Default is \`TRUE\`.

## Value

Data frame with cells as rows and features as columns.

## Author

Dongqiang Zeng

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("Seurat", quietly = TRUE)) {
  pbmc <- SeuratObject::pbmc_small
  vars <- c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3")
  eset <- extract_sc_data(sce = pbmc, vars = vars, assay = "RNA")
}
} # }
```
