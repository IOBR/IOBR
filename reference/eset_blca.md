# TCGA-BLCA Bladder Cancer Expression Data

Gene expression data from bladder cancer patients in TCGA-BLCA. Contains
raw count matrix suitable for differential expression analysis.

## Usage

``` r
data(eset_blca)
```

## Format

A numeric matrix with genes as rows and samples as columns:

- Rows:

  Ensembl gene identifiers (e.g., ENSG00000000003)

- Columns:

  TCGA sample identifiers (e.g., TCGA-2F-A9KO)

- Values:

  Raw expression counts (integer values)

## Source

The Cancer Genome Atlas Bladder Urothelial Carcinoma (TCGA-BLCA)

## References

Cancer Genome Atlas Research Network. Comprehensive molecular
characterization of urothelial bladder carcinoma. Nature 507, 315-322
(2014). doi:10.1038/nature12965
