# Use quanTIseq to Deconvolute a Gene Expression Matrix

Deconvolutes gene expression data to estimate immune cell fractions
using the quanTIseq method. Source code from
https://github.com/FFinotello/quanTIseq.

## Usage

``` r
deconvolute_quantiseq.default(
  mix.mat,
  arrays = FALSE,
  signame = "TIL10",
  tumor = FALSE,
  mRNAscale = TRUE,
  method = c("lsei", "hampel", "huber", "bisquare"),
  rmgenes = "unassigned"
)
```

## Arguments

- mix.mat:

  Data frame or matrix. Gene expression matrix with gene symbols on the
  first column and sample IDs on the first row. Expression data must be
  on non-log scale (TPM for RNA-seq or expression values for
  microarrays).

- arrays:

  Logical. Whether expression data are from microarrays. Default is
  FALSE. If TRUE, the rmgenes parameter is set to "none".

- signame:

  Character. Name of the signature matrix. Currently only "TIL10" is
  available. Default is "TIL10".

- tumor:

  Logical. Whether expression data are from tumor samples. If TRUE,
  signature genes with high expression in tumor samples are removed.
  Default is FALSE.

- mRNAscale:

  Logical. Whether cell fractions must be scaled to account for
  cell-type-specific mRNA content. Default is TRUE.

- method:

  Character. Deconvolution method: "hampel", "huber", "bisquare" for
  robust regression, or "lsei" for constrained least squares. Default is
  "lsei".

- rmgenes:

  Character. Genes to remove: "unassigned" (default), "default", "none",
  or "path".

## Value

Data frame with cell fractions for each sample.

## References

F. Finotello, C. Mayer, C. Plattner, G. Laschober, D. Rieder, H. Hackl,
A. Krogsdam, W. Posch, D. Wilflingseder, S. Sopper, M. Jsselsteijn, D.
Johnsons, Y. Xu, Y. Wang, M. E. Sanders, M. V. Estrada, P.
Ericsson-Gonzalez, J. Balko, N. F. de Miranda, Z. Trajanoski.
"quanTIseq: quantifying immune contexture of human tumors". bioRxiv
223180. https://doi.org/10.1101/223180.

## Author

Finotello F, et al. (adapted for IOBR)

## Examples

``` r
quantiseq_data <- load_data("quantiseq_data")
#> ℹ Loading cached data: "quantiseq_data"
if (!is.null(quantiseq_data)) {
  common_genes <- rownames(quantiseq_data$TIL10_signature)
  tpm_matrix <- as.data.frame(matrix(
    abs(rnorm(length(common_genes) * 5, mean = 5, sd = 2)),
    nrow = length(common_genes), ncol = 5
  ))
  rownames(tpm_matrix) <- common_genes
  colnames(tpm_matrix) <- paste0("Sample", 1:5)
  results <- deconvolute_quantiseq.default(mix.mat = tpm_matrix)
  if (!is.null(results)) head(results)
}
#> ℹ Running quanTIseq deconvolution module
#> ℹ Loading cached data: "quantiseq_data"
#> ℹ Gene expression normalization and re-annotation (arrays: FALSE)
#> ℹ Loading cached data: "quantiseq_data"
#> ℹ Removing 17 noisy genes
#> ℹ Signature genes found in data set: 153/153 (100%)
#> ℹ Mixture deconvolution (method: lsei)
#> ✔ Deconvolution successful!
#>          Sample B.cells Macrophages.M1 Macrophages.M2 Monocytes Neutrophils
#> Sample1 Sample1       0              0              0         0           0
#> Sample2 Sample2       1              0              0         0           0
#> Sample3 Sample3       0              0              0         0           0
#> Sample4 Sample4       0              1              0         0           0
#> Sample5 Sample5       0              0              1         0           0
#>         NK.cells T.cells.CD4 T.cells.CD8 Tregs Dendritic.cells        Other
#> Sample1        1           0           0     0               0 0.000000e+00
#> Sample2        0           0           0     0               0 0.000000e+00
#> Sample3        1           0           0     0               0 0.000000e+00
#> Sample4        0           0           0     0               0 1.376677e-14
#> Sample5        0           0           0     0               0 0.000000e+00
```
