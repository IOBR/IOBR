# GSE62254 Gastric Cancer Expression Data

Gene expression data from gastric cancer patients in the GSE62254
dataset. Contains raw count matrix from RNA-seq or microarray
experiments.

## Usage

``` r
data(eset_gse62254)
```

## Format

A numeric matrix with genes as rows and samples as columns:

- Rows:

  Ensembl gene identifiers (e.g., ENSG00000000003)

- Columns:

  Sample identifiers (e.g., TCGA-2F-A9KQ)

- Values:

  Expression counts or intensities

## Source

GEO accession GSE62254

## References

Cristescu R et al. (2015) Molecular analysis of gastric cancer
identifies subtypes associated with distinct clinical outcomes. Nat Med
21:449-456. doi:10.1038/nm.3850

## Examples

``` r
data(eset_gse62254)
head(eset_gse62254)
#>           GSM1523727 GSM1523728 GSM1523729
#> 1007_s_at  3.2176645  3.0624323  3.0279131
#> 1053_at    2.4050109  2.4394879  2.2442708
#> 117_at     1.4933412  1.8067380  1.5959665
#> 121_at     2.1965561  2.2812181  2.1865556
#> 1255_g_at  0.8698382  0.9502466  0.8125414
#> 1294_at    2.1203222  2.1587338  2.2069812
```
