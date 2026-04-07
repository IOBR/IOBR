# Generate Reference Signature Matrix Using DESeq2

Uses DESeq2 to perform differential expression analysis across cell
types, identifies significantly expressed genes, and creates a reference
signature matrix from median expression levels.

## Usage

``` r
generateRef_DEseq2(dds, pheno, FDR = 0.05, dat)
```

## Arguments

- dds:

  Matrix. Raw count data from RNA-seq.

- pheno:

  Character vector. Cell type classes for samples.

- FDR:

  Numeric. Threshold for adjusted p-values. Default is 0.05.

- dat:

  Matrix. Normalized expression data (e.g., FPKM, TPM) for calculating
  median expression.

## Value

List containing: - \`reference_matrix\`: Data frame of median expression
for significant genes across cell types. - \`G\`: Optimal number of
probes minimizing condition number. - \`condition_number\`: Minimum
condition number. - \`whole_matrix\`: Full median expression matrix.

## Examples

``` r
set.seed(123)
dds <- matrix(sample(0:1000, 2000, replace = TRUE), nrow = 100, ncol = 20)
colnames(dds) <- paste("Sample", 1:20, sep = "_")
rownames(dds) <- paste("Gene", 1:100, sep = "_")
pheno <- rep(c("Type1", "Type2"), each = 10)
dat <- matrix(runif(2000), nrow = 100, ncol = 20)
rownames(dat) <- rownames(dds)
colnames(dat) <- colnames(dds)
# \donttest{
result <- generateRef_DEseq2(dds = dds, pheno = pheno, FDR = 0.05, dat = dat)
#> converting counts to integer mode
#> Warning: some variables in design formula are characters, converting to factors
#> using pre-existing size factors
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
#> using pre-existing size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
#> final dispersion estimates
#> fitting model and testing
print(result$reference_matrix)
#> [1] Type1 Type2
#> <0 rows> (or 0-length row.names)
# }
```
