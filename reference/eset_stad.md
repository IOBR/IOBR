# Toy STAD expression matrix

A subset of TCGA-STAD RNA-seq count data used in the IOBR vignette and
unit tests. Rows correspond to Ensembl gene IDs and columns correspond
to TCGA sample barcodes. Values are raw (unnormalized) read counts.

## Usage

``` r
data(eset_stad)
```

## Format

A numeric matrix with genes in rows and samples in columns.
