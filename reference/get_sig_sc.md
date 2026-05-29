# Extract Top Marker Genes from Single-Cell Differential Results

Selects the top N marker genes per cluster from a ranked differential
expression result table.

## Usage

``` r
get_sig_sc(
  deg,
  cluster = "cluster",
  gene = "gene",
  avg_log2FC = "avg_log2FC",
  n = 100
)
```

## Arguments

- deg:

  Data frame or matrix. Ranked marker statistics.

- cluster:

  Character. Column name containing cluster identifiers. Default is
  \`"cluster"\`.

- gene:

  Character. Column name containing gene identifiers. Default is
  \`"gene"\`.

- avg_log2FC:

  Character. Column name for average log2 fold change. Default is
  \`"avg_log2FC"\`.

- n:

  Integer. Number of top markers per cluster. Default is 100.

## Value

List of character vectors; each element contains the top N genes for a
cluster.

## Examples

``` r
# Simulate marker data
set.seed(123)
sim_deg <- data.frame(
  cluster = rep(c("A", "B"), each = 50),
  gene = paste0("Gene", 1:100),
  avg_log2FC = rnorm(100, 2, 1)
)

# Extract top 5 markers per cluster
markers <- get_sig_sc(sim_deg, n = 5)
print(markers)
#> $A
#> [1] "Gene3"  "Gene6"  "Gene16" "Gene30" "Gene44"
#> 
#> $B
#> [1] "Gene54" "Gene56" "Gene70" "Gene97" "Gene98"
#> 
```
