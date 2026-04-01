# GRCh38 Human Genome Annotation

Gene annotation for human GRCh38/hg38 genome assembly, used for
expression normalization and gene identifier mapping in IOBR functions.

## Usage

``` r
data(anno_grch38)
```

## Format

A data frame with 60668 rows and 11 columns:

- id:

  Ensembl gene identifier.

- eff_length:

  Effective gene length for TPM calculation.

- gc:

  GC content proportion.

- entrez:

  Entrez Gene ID.

- symbol:

  Official gene symbol.

- chr:

  Chromosome name.

- start:

  Genomic start position.

- end:

  Genomic end position.

- strand:

  Genomic strand.

- biotype:

  Gene biotype classification.

- description:

  Gene description.

## Source

Ensembl release 104 (GRCh38.p13)
