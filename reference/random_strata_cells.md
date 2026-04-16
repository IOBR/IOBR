# Stratified Random Sampling of Cells

Performs stratified random sampling of cells from single-cell data,
ensuring proportional representation of each cell type while respecting
minimum and maximum count constraints.

## Usage

``` r
random_strata_cells(
  input,
  group,
  proportion = 0.1,
  minimum_count_include = 300,
  minimum_count = 200,
  maximum_count = 1000,
  sub_cluster = NULL,
  cell_type = NULL
)
```

## Arguments

- input:

  A data frame or Seurat object containing cell annotations.

- group:

  Character string specifying the column name for cell type grouping.

- proportion:

  Numeric value between 0 and 1 specifying the sampling proportion.
  Default is 0.1.

- minimum_count_include:

  Integer specifying the minimum count threshold for a cell type to be
  included in sampling. Default is 300.

- minimum_count:

  Integer specifying the minimum number of cells to sample per cell
  type. Default is 200.

- maximum_count:

  Integer specifying the maximum number of cells to sample per cell
  type. Default is 1000.

- sub_cluster:

  Optional character string specifying a sub-cluster column for
  filtering. Default is NULL.

- cell_type:

  Optional character string specifying the cell type value to filter
  when \`sub_cluster\` is provided. Default is NULL.

## Value

A data frame containing the sampled cells with preserved cell type
proportions.

## Examples

``` r
# \donttest{
# Sample cells from a data frame
sampled_cells <- random_strata_cells(
  input = cell_annotations,
  group = "cell_type",
  proportion = 0.1,
  minimum_count_include = 300,
  minimum_count = 200,
  maximum_count = 1000
)
#> Error: object 'cell_annotations' not found

# Sample cells from a Seurat object
sampled_cells <- random_strata_cells(
  input = seurat_object,
  group = "seurat_clusters",
  proportion = 0.2
)
#> Error: object 'seurat_object' not found
# }
```
