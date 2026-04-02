# Generate Reference Gene Matrix from RNA-seq DEGs

Uses DESeq2 to identify differentially expressed genes and create a
reference matrix from median expression levels across cell types.

## Usage

``` r
generateRef_rnaseq(dds, pheno, mode = "oneVSothers", FDR = 0.05, dat)
```

## Arguments

- dds:

  Matrix. Raw count data from RNA-seq.

- pheno:

  Character vector. Cell type classes.

- mode:

  Character. DEG identification mode: \`"oneVSothers"\` or \`"pairs"\`.
  Default is \`"oneVSothers"\`.

- FDR:

  Numeric. Threshold for adjusted p-values. Default is 0.05.

- dat:

  Matrix. Normalized expression data (e.g., FPKM, TPM).

## Value

List containing: - \`reference_matrix\`: Data frame of median expression
for significant genes across cell types. - \`G\`: Optimal number of
probes minimizing condition number. - \`condition_number\`: Minimum
condition number. - \`whole_matrix\`: Full median expression matrix.

## Author

Rongfang Shen

## Examples

``` r
dds <- matrix(rpois(200 * 10, lambda = 10), ncol = 10)
rownames(dds) <- paste("Gene", 1:200, sep = "_")
colnames(dds) <- paste("Sample", 1:10, sep = "_")
pheno <- sample(c("Type1", "Type2", "Type3"), 10, replace = TRUE)
dat <- matrix(rnorm(200 * 10), ncol = 10)
rownames(dat) <- rownames(dds)
colnames(dat) <- colnames(dds)
# \donttest{
results <- generateRef_rnaseq(dds = dds, pheno = pheno, FDR = 0.05, dat = dat)
#> converting counts to integer mode
#> Warning: some variables in design formula are characters, converting to factors
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> converting counts to integer mode
#> Warning: some variables in design formula are characters, converting to factors
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
#> final dispersion estimates
#> fitting model and testing
#> converting counts to integer mode
#> Warning: some variables in design formula are characters, converting to factors
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
#> final dispersion estimates
#> fitting model and testing
# }
```
