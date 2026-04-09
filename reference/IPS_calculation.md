# Calculate Immunophenoscore (IPS)

Calculates Immunophenoscore (IPS) from gene expression data. IPS is a
composite score measuring immunophenotype based on four major
categories: MHC molecules, immunomodulators, effector cells, and
suppressor cells.

## Usage

``` r
IPS_calculation(project = NULL, eset, plot = FALSE)
```

## Arguments

- project:

  Character string for project identifier. Default is NULL.

- eset:

  Gene expression matrix with official human gene symbols (HGNC) as
  rownames. Expression values should be log2(TPM+1) or will be
  transformed if max value \> 100.

- plot:

  Logical. Whether to generate immunophenogram plots. Default is FALSE.

## Value

Data frame containing:

- MHC:

  MHC molecules score

- EC:

  Effector cells score

- SC:

  Suppressor cells score

- CP:

  Checkpoints/Immunomodulators score

- AZ:

  Aggregate score (sum of MHC, CP, EC, SC)

- IPS:

  Immunophenoscore (0-10 scale)

## Examples

``` r
# IPS requires gene symbols as rownames
# Create a simple example with gene symbols
example_genes <- c(
  "HLA-A", "HLA-B", "HLA-C", "CD274", "PDCD1", "CTLA4",
  "CD8A", "CD8B", "GZMB", "PRF1", "FOXP3", "IL10"
)
sim_eset <- as.data.frame(matrix(
  rnorm(length(example_genes) * 10, mean = 5, sd = 2),
  nrow = length(example_genes), ncol = 10
))
rownames(sim_eset) <- example_genes
colnames(sim_eset) <- paste0("Sample", 1:10)
ips_result <- IPS_calculation(eset = sim_eset, project = "Example", plot = FALSE)
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "ips_gene_set"
head(ips_result)
#>         ProjectID         MHC  EC  SC  CP          AZ IPS
#> Sample1   Example  0.84758794 NaN NaN NaN  0.84758794   3
#> Sample2   Example  0.17473579 NaN NaN NaN  0.17473579   1
#> Sample3   Example  0.29079710 NaN NaN NaN  0.29079710   1
#> Sample4   Example  0.19513209 NaN NaN NaN  0.19513209   1
#> Sample5   Example -0.08050166 NaN NaN NaN -0.08050166   0
#> Sample6   Example -0.20960706 NaN NaN NaN -0.20960706   0
```
