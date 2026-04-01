# TCGA-STAD Tumor Microenvironment Signature Scores

Pre-calculated tumor microenvironment (TME) signature scores for TCGA
stomach adenocarcinoma (STAD) samples. Contains expression levels of key
TME-related genes and potentially immune/stomal signature scores.

## Usage

``` r
data(eset_tme_stad)
```

## Format

A numeric matrix with genes/signatures as rows and samples as columns:

- Rows:

  Gene symbols or signature names (e.g., B2M, HLA-B, HSPB1)

- Columns:

  TCGA sample identifiers (e.g., TCGA-3M-AB46-01A)

- Values:

  Normalized expression values (appears to be log2-scale)

## Source

The Cancer Genome Atlas Stomach Adenocarcinoma (TCGA-STAD)

## References

Cancer Genome Atlas Research Network. Comprehensive molecular
characterization of gastric adenocarcinoma. Nature. 2014;513:202-209.
doi:10.1038/nature13480
