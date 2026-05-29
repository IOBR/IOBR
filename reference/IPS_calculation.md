# Calculate Immunophenoscore (IPS)

Calculates Immunophenoscore (IPS) from gene expression data. IPS is a
composite score measuring immunophenotype based on four major
categories: MHC molecules, immunomodulators, effector cells, and
suppressor cells.

## Usage

``` r
IPS_calculation(eset, project = NULL, plot = FALSE)
```

## Arguments

- eset:

  Gene expression matrix with official human gene symbols (HGNC) as
  rownames. Expression values should be log2(TPM+1) or will be
  transformed if max value \> 100.

- project:

  Character string for project identifier. Default is NULL.

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
if (interactive()) {
  example_genes <- c(
    "HLA-A", "HLA-B", "HLA-C", "CD274", "PDCD1", "CTLA4",
    "CD8A", "CD8B", "GZMB", "PRF1", "FOXP3", "IL10"
  )
  sim_eset <- as.data.frame(matrix(
    rnorm(length(example_genes) * 5, mean = 5, sd = 2),
    nrow = length(example_genes), ncol = 5
  ))
  rownames(sim_eset) <- example_genes
  colnames(sim_eset) <- paste0("Sample", 1:5)
  ips_result <- IPS_calculation(eset = sim_eset, project = "Example", plot = FALSE)
  if (!is.null(ips_result)) head(ips_result)
}
```
