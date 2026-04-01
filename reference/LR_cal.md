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
# \donttest{
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
if (requireNamespace("easier", quietly = TRUE)) {
  lr <- LR_cal(eset = sim_eset, data_type = "tpm")
  head(lr)
}
#> 
#> 
#> 
#> LR signature genes found in data set: 11/644 (1.7%)
#> Ligand-Receptor pair weights computed 
#>        ID CXCL12_CXCR4   IL6_IL6R TNF_TNFRSF1A_TNFRSF21_TRAF2
#> 1 Sample1     2.368152  2.4361612                    2.341010
#> 2 Sample2     2.759510 -0.2755479                    1.052273
#> 3 Sample3     2.057308  2.7174774                    2.716119
#> 4 Sample4     2.426236  2.1919763                    2.178430
#> 5 Sample5     2.328362  2.6795345                    1.777969
#> 6 Sample6     2.582983  2.5136527                    2.255756
# }
```
