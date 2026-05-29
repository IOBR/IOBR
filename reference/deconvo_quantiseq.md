# Deconvolve Using quanTIseq

quanTIseq deconvolution for RNA-seq immune cell fractions.

## Usage

``` r
deconvo_quantiseq(eset, project = NULL, tumor, arrays, scale_mrna)
```

## Arguments

- eset:

  Gene expression matrix.

- project:

  Optional project name. Default is \`NULL\`.

- tumor:

  Logical: tumor samples. Must be specified.

- arrays:

  Logical: microarray data. Must be specified.

- scale_mrna:

  Logical: correct for mRNA content. Must be specified.

## Value

Data frame with quanTIseq cell fractions. Columns suffixed with
\`\_quantiseq\`.

## Author

Dongqiang Zeng

## Examples

``` r
if (interactive()) {
  quantiseq_data <- load_data("quantiseq_data")
  if (!is.null(quantiseq_data)) {
    set.seed(123)
    n_sig <- nrow(quantiseq_data$TIL10_signature)
    sim_eset <- matrix(rnorm(n_sig * 2), n_sig, 2)
    rownames(sim_eset) <- rownames(quantiseq_data$TIL10_signature)
    colnames(sim_eset) <- paste0("Sample", 1:2)
    result <- deconvo_quantiseq(eset = sim_eset, project = "Example", tumor = TRUE,
                                arrays = FALSE, scale_mrna = FALSE)
    if (!is.null(result)) head(result)
  }
}
```
