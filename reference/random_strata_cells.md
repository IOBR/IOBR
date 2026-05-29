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
# Create simulated cell annotation data
set.seed(123)
sim_cells <- data.frame(
  cell_id = paste0("Cell", 1:500),
  cell_type = sample(c("T_cell", "B_cell", "NK_cell", "Macrophage"), 500, replace = TRUE)
)
# Sample cells with stratified random sampling
sampled <- random_strata_cells(
  input = sim_cells,
  group = "cell_type",
  proportion = 0.2,
  minimum_count_include = 50,
  minimum_count = 20,
  maximum_count = 100
)
#> ℹ Cell type counts before sampling:
#> ℹ Cell types included in sampling:
#> ℹ Initial sample sizes (proportion = 0.2):
#> ℹ Cell type counts after sampling:
if (!is.null(sampled)) head(sampled)
#>     cell_id cell_type
#> 48   Cell48    B_cell
#> 59   Cell59    B_cell
#> 81   Cell81    B_cell
#> 92   Cell92    B_cell
#> 103 Cell103    B_cell
#> 120 Cell120    B_cell
```
