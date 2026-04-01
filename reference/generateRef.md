# Generate Reference Signature Matrix

Generates a reference signature matrix for cell types based on
differential expression analysis. Supports both limma for normalized
data and DESeq2 for raw count data.

## Usage

``` r
generateRef(dds, pheno, FDR = 0.05, dat, method = "limma")
```

## Arguments

- dds:

  Matrix. Raw count data from RNA-seq. Required if \`method =
  "DESeq2"\`.

- pheno:

  Character vector. Cell type class of the samples.

- FDR:

  Numeric. Genes with BH adjusted p-value \< FDR are considered
  significant. Default is 0.05.

- dat:

  Matrix or data frame. Normalized transcript quantification data (e.g.,
  FPKM, TPM).

- method:

  Character. Method for differential expression: \`"limma"\` or
  \`"DESeq2"\`. Default is \`"limma"\`.

## Value

List containing: - \`reference_matrix\`: Data frame of median expression
for significant genes across cell types. - \`G\`: Optimal number of
probes minimizing condition number. - \`condition_number\`: Minimum
condition number. - \`whole_matrix\`: Full median expression matrix.

## Examples

``` r
expressionData <- matrix(runif(1000 * 4, min = 0, max = 10), ncol = 4)
rownames(expressionData) <- paste("Gene", 1:1000, sep = "_")
colnames(expressionData) <- paste("Sample", 1:4, sep = "_")

phenotype <- c("celltype1", "celltype2", "celltype1", "celltype2")

rawCountData <- matrix(sample(1:100, 1000 * 4, replace = TRUE), ncol = 4)
rownames(rawCountData) <- paste("Gene", 1:1000, sep = "_")
colnames(rawCountData) <- paste("Sample", 1:4, sep = "_")

# \donttest{
result <- generateRef(
  dds = rawCountData, pheno = phenotype,
  FDR = 0.05, dat = expressionData, method = "DESeq2"
)
#> ℹ Running differentially expressed genes using DESeq2
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
