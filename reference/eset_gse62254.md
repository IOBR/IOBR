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
