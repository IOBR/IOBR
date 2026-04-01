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
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
eset <- eset[1:500, 1:3]
res <- deconvo_timer(eset = eset, project = "stad", indications = rep("stad", ncol(eset)))
} # }
```
