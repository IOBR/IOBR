# Deconvolve Using TIMER

TIMER deconvolution for cancer-specific immune estimation.

## Usage

``` r
deconvo_timer(eset, project = NULL, indications = NULL)
```

## Arguments

- eset:

  Gene expression matrix.

- project:

  Optional project name. Default is \`NULL\`.

- indications:

  Cancer type for each sample (e.g., \`"brca"\`, \`"stad"\`). Must match
  number of columns in \`eset\`.

## Value

Data frame with TIMER cell fractions. Columns suffixed with \`\_TIMER\`.

## Author

Dongqiang Zeng

## Examples

``` r
if (FALSE) { # \dontrun{
immune <- load_data("immuneCuratedData")
cancer_genes <- load_data("cancer_type_genes")
if (!is.null(immune) && !is.null(cancer_genes)) {
  set.seed(123)
  genes <- unique(c(head(rownames(immune$genes), 100), cancer_genes[["stad"]]))
  sim_eset <- matrix(rnorm(length(genes) * 2), length(genes), 2)
  rownames(sim_eset) <- genes
  colnames(sim_eset) <- paste0("Sample", 1:2)
  result <- deconvo_timer(eset = sim_eset, project = "TCGA-STAD",
                          indications = rep("stad", 2))
  if (!is.null(result)) head(result)
}
} # }
```
