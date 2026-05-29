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
# Simulate data
set.seed(123)
immune <- load_data("immuneCuratedData")
#> ℹ Loading cached data: "immuneCuratedData"
cancer_genes <- load_data("cancer_type_genes")
#> ℹ Loading cached data: "cancer_type_genes"
if (!is.null(immune) && !is.null(cancer_genes)) {
  gene_names <- unique(c(head(rownames(immune$genes), 50), cancer_genes[["stad"]]))
  sample_names <- paste0("Sample", 1:5)
  n_genes <- length(gene_names)
  expr <- matrix(runif(n_genes * 5, 1, 100), n_genes, 5,
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
#> ℹ Enter batch mode
#> ℹ Loading immune gene expression
#> ℹ Loading cached data: "immuneCuratedData"
#> ℹ Outlier genes: A2BP1 AADACL2 AASDHPPT ABCB8 AIM2 BTLA C10orf54 CD247 GPR18 GZMH HCLS1 HS3ST3B1 IL15 ITGAM ITK LCP2 LILRB1 MIR155HG PLA2G7 S100A12 STAT5A TIFAB TREM2 TSC22D3
#> ℹ Removing batch effects for stad
#> ℹ Loading cached data: "cancer_type_genes"
#>               Sample1    Sample2    Sample3    Sample4    Sample5
#> B_cell     0.09675834 0.11148673 0.09241700 0.10200825 0.09059373
#> T_cell.CD4 0.12579438 0.12470773 0.12391277 0.11720438 0.12227773
#> T_cell.CD8 0.19159022 0.18977545 0.19816059 0.19191750 0.19316088
#> Neutrophil 0.10771851 0.10756189 0.10473387 0.10691379 0.10657790
#> Macrophage 0.05270592 0.03984843 0.04274168 0.04611653 0.03743177
#> DC         0.46616920 0.47245314 0.47769228 0.47863637 0.48888020
```
