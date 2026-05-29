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
if (interactive()) {
  ips_genes <- load_data("ips_gene_set")
  if (!is.null(ips_genes)) {
    set.seed(123)
    sim_eset <- matrix(rnorm(nrow(ips_genes) * 2), nrow(ips_genes), 2)
    rownames(sim_eset) <- ips_genes$GENE
    colnames(sim_eset) <- paste0("Sample", 1:2)
    result <- deconvo_ips(eset = sim_eset, project = "Example")
    if (!is.null(result)) head(result)
  }
}
```
