# Calculate Ligand-Receptor Interaction Scores

Quantifies ligand-receptor interactions in the tumor microenvironment
from bulk gene expression data using the easier package. This function
processes raw counts or TPM data and computes interaction scores for
each sample.

## Usage

``` r
LR_cal(
  eset,
  data_type = c("count", "tpm"),
  id_type = "ensembl",
  cancer_type = "pancan"
)
```

## Arguments

- eset:

  Gene expression matrix with genes as rows and samples as columns.

- data_type:

  Type of input data. Options are \`"count"\` or \`"tpm"\`. If
  \`"count"\`, data will be converted to TPM before analysis.

- id_type:

  Type of gene identifier. Default is \`"ensembl"\`.

- cancer_type:

  Character string specifying the cancer type for easier. Default is
  \`"pancan"\` for pan-cancer analysis.

## Value

Data frame containing ligand-receptor interaction scores with sample IDs
as row names.

## References

Lapuente-Santana, van Genderen, M., Hilbers, P., Finotello, F., &
Eduati, F. (2021). Interpretable systems biomarkers predict response to
immune-checkpoint inhibitors. Patterns (New York, N.Y.), 2(8), 100293.
https://doi.org/10.1016/j.patter.2021.100293

## Examples

``` r
# LR_cal requires HGNC gene symbols as rownames
# Create a simple example with gene symbols
example_genes <- c(
  "TGFB1", "EGFR", "VEGFA", "PDGFB", "FGF2", "CXCL12",
  "CXCR4", "IL6", "IL6R", "TNF", "TNFRSF1A", "IFNG"
)
sim_eset <- as.data.frame(matrix(
  rnorm(length(example_genes) * 10, mean = 5, sd = 2),
  nrow = length(example_genes), ncol = 10
))
rownames(sim_eset) <- example_genes
colnames(sim_eset) <- paste0("Sample", 1:10)
if (FALSE) { # \dontrun{
if (requireNamespace("easier", quietly = TRUE)) {
  tryCatch({
    lr <- LR_cal(eset = sim_eset, data_type = "tpm")
    head(lr)
  }, error = function(e) {
    message("Example skipped: could not download ExperimentHub data")
  })
}
} # }
```
