# Grouped gene signatures for IOBR analysis

A named list that organizes gene signatures into functional or
biological categories. Each element of the list is a character vector
containing the names of gene signatures defined in
`signature_collection`. A total of 43 signature groups are included,
covering tumour intrinsic pathways, immune-related processes, stromal
activity, TME characteristics and immuno-oncology biomarkers. These
groups are used in IOBR to conveniently select sets of signatures for
scoring and visualization.

## Usage

``` r
data(sig_group)
```

## Format

A named list of length 43. Each element is a character vector of
signature names. Representative groups include:

- tumor_signature:

  Signatures related to intrinsic tumour biology such as cell cycle, DNA
  damage repair and histone regulation.

- EMT:

  Epithelial–mesenchymal transition (EMT)–associated signatures.

- io_biomarkers:

  Immuno-oncology biomarker–related signatures.

- immu_microenvironment:

  Immune microenvironment–related signatures.

- immu_suppression:

  Immune suppression–related signatures.

- immu_exclusion:

  Signatures associated with immune exclusion and stromal barriers.

- TCR_BCR:

  T-cell and B-cell receptor pathway signatures.

- tme_signatures1:

  Tumour microenvironment signature panel (set 1).

- tme_signatures2:

  Tumour microenvironment signature panel (set 2).

- ...:

  Additional groups are included (43 total), but not listed individually
  here; all groups follow the same structure.
