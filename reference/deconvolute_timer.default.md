# Deconvolute Tumor Microenvironment Using TIMER

Performs deconvolution of the tumor microenvironment using the TIMER
algorithm. Processes multiple cancer datasets, removes batch effects,
and estimates immune cell type abundances.

## Usage

``` r
deconvolute_timer.default(args)
```

## Arguments

- args:

  List or environment containing parameters:

  outdir

  :   Character. Output directory path.

  batch

  :   Character. File containing paths to expression data and cancer
      types.

## Value

Matrix of abundance scores for different immune cell types across
multiple cancer samples.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
immune <- load_data("immuneCuratedData")
cancer_genes <- load_data("cancer_type_genes")
if (!is.null(immune) && !is.null(cancer_genes)) {
  gene_names <- unique(c(head(rownames(immune$genes), 30), cancer_genes[["stad"]]))
  sample_names <- paste0("Sample", 1:2)
  n_genes <- length(gene_names)
  expr <- matrix(runif(n_genes * 2, 1, 100), n_genes, 2,
                 dimnames = list(gene_names, sample_names))
  tf <- tempfile(fileext = ".tsv")
  write.table(as.data.frame(expr), tf, sep = "\t", quote = FALSE)
  batch_tf <- tempfile(fileext = ".csv")
  write.table(data.frame(path = tf, type = "stad"), batch_tf,
              sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
  args <- list(outdir = tempdir(), batch = batch_tf)
  result <- deconvolute_timer.default(args)
  if (!is.null(result)) head(result)
}
} # }
```
