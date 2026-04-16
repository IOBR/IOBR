# Generate Reference Matrix from Seurat Object

Generates reference gene expression data from a Seurat object by
identifying marker genes for each cell type and aggregating expression
data.

## Usage

``` r
generateRef_seurat(
  sce,
  celltype = NULL,
  proportion = NULL,
  assay_deg = "RNA",
  slot_deg = "data",
  adjust_assay = FALSE,
  assay_out = "RNA",
  slot_out = "data",
  verbose = FALSE,
  only.pos = TRUE,
  n_ref_genes = 50,
  logfc.threshold = 0.15,
  test.use = "wilcox"
)
```

## Arguments

- sce:

  Seurat object containing single-cell RNA-seq data.

- celltype:

  Character. Cell type column name in metadata. Default is \`NULL\`
  (uses default identity).

- proportion:

  Numeric. Proportion of cells to randomly select for analysis. Default
  is \`NULL\` (use all cells).

- assay_deg:

  Character. Assay for finding markers. Default is \`"RNA"\`.

- slot_deg:

  Character. Slot for finding markers. Default is \`"data"\`.

- adjust_assay:

  Logical. Whether to adjust assay for SCT. Default is \`FALSE\`.

- assay_out:

  Character. Assay for output. Default is \`"RNA"\`.

- slot_out:

  Character. Slot for output. Default is \`"data"\`.

- verbose:

  Logical. Print verbose messages. Default is \`FALSE\`.

- only.pos:

  Logical. Return only positive markers. Default is \`TRUE\`.

- n_ref_genes:

  Integer. Number of reference genes per cell type. Default is 50.

- logfc.threshold:

  Numeric. Log fold change threshold. Default is 0.15.

- test.use:

  Character. Statistical test for marker identification. Default is
  \`"wilcox"\`.

## Value

Matrix containing aggregated expression data for reference genes.

## Author

Dongqiang Zeng

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("Seurat", quietly = TRUE)) {
  # Requires a Seurat object with sufficient cells and markers
  sm <- generateRef_seurat(sce = seurat_obj, celltype = "cell_type", slot_out = "data")
}
} # }
```
