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
if (FALSE) { # \dontrun{
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
eset <- eset[1:500, 1:3]
res <- deconvo_quantiseq(
  eset = eset, project = "stad", tumor = TRUE,
  arrays = FALSE, scale_mrna = FALSE
)
} # }
```
