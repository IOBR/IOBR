# Calculate ESTIMATE Scores

Calculates stromal, immune, and ESTIMATE scores from gene expression.

## Usage

``` r
deconvo_estimate(eset, project = NULL, platform = "affymetrix")
```

## Arguments

- eset:

  Gene expression matrix with gene symbols.

- project:

  Optional project name. Default is \`NULL\`.

- platform:

  Platform type: \`"affymetrix"\` or \`"illumina"\`. Default is
  \`"affymetrix"\`.

## Value

Data frame with ESTIMATE scores. Columns suffixed with \`\_estimate\`.

## Author

Dongqiang Zeng

## Examples

``` r
if (interactive()) {
  common_genes <- load_data("common_genes")
  if (!is.null(common_genes)) {
    set.seed(123)
    sim_eset <- matrix(rnorm(nrow(common_genes) * 2), nrow(common_genes), 2)
    rownames(sim_eset) <- common_genes$GeneSymbol
    colnames(sim_eset) <- paste0("Sample", 1:2)
    result <- deconvo_estimate(sim_eset, project = "TCGA-STAD")
    if (!is.null(result)) head(result)
  }
}
```
