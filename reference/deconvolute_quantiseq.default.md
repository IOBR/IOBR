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
# \donttest{
lm22 <- load_data("lm22")
common_genes <- rownames(lm22)[1:500]
tpm_matrix <- as.data.frame(matrix(
  rnorm(length(common_genes) * 5, mean = 5, sd = 2),
  nrow = length(common_genes), ncol = 5
))
rownames(tpm_matrix) <- common_genes
colnames(tpm_matrix) <- paste0("Sample", 1:5)
results <- deconvolute_quantiseq.default(mix.mat = tpm_matrix)
#> ℹ Running quanTIseq deconvolution module
#> ℹ Gene expression normalization and re-annotation (arrays: FALSE)
#> ℹ Removing 17 noisy genes
#> ℹ Signature genes found in data set: 41/153 (26.8%)
#> ℹ Mixture deconvolution (method: lsei)
#> ✔ Deconvolution successful!
# }
```
