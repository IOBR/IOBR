# Calculate Immunophenoscore (IPS)

Calculates immune phenotype scores from gene expression data.

## Usage

``` r
deconvo_ips(eset, project = NULL, plot = FALSE)
```

## Arguments

- eset:

  Gene expression matrix.

- project:

  Optional project name. Default is \`NULL\`.

- plot:

  Logical: generate visualization. Default is \`FALSE\`.

## Value

Data frame with IPS scores. Columns suffixed with \`\_IPS\`.

## Author

Dongqiang Zeng

## Examples

``` r
ips_genes <- load_data("ips_gene_set")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "ips_gene_set"
if (!is.null(ips_genes)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(nrow(ips_genes) * 3), nrow(ips_genes), 3)
  rownames(sim_eset) <- ips_genes$GENE
  colnames(sim_eset) <- paste0("Sample", 1:3)
  
  # Run calculation
  result <- deconvo_ips(eset = sim_eset, project = "Example")
  if (!is.null(result)) head(result)
}
#> ℹ Running IPS calculation
#> ℹ Loading cached data: "ips_gene_set"
#>        ID ProjectID ProjectID_IPS    MHC_IPS      EC_IPS      SC_IPS     CP_IPS
#> 1 Sample1   Example       Example 0.09554254 -0.02191850  0.07525168 -0.5033694
#> 2 Sample2   Example       Example 0.20070607 -0.06235335 -0.14211063 -0.3425485
#> 3 Sample3   Example       Example 0.23776174 -0.09747729 -0.20235106  0.4555573
#>       AZ_IPS IPS_IPS
#> 1 -0.3544937       0
#> 2 -0.3463064       0
#> 3  0.3934907       1
```
