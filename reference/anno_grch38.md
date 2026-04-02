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

## Examples

``` r
data(anno_grch38)
head(anno_grch38)
#>                id eff_length        gc entrez   symbol chr     start       end
#> 1 ENSG00000000003       4536 0.3992504   7105   TSPAN6   X 100627109 100639991
#> 2 ENSG00000000005       1476 0.4241192  64102     TNMD   X 100584802 100599885
#> 3 ENSG00000000419       9276 0.4252911   8813     DPM1  20  50934867  50958555
#> 4 ENSG00000000457       6883 0.4117391  57147    SCYL3   1 169849631 169894267
#> 5 ENSG00000000460       5970 0.4298157  55732 C1orf112   1 169662007 169854080
#> 6 ENSG00000000938       3382 0.5644589   2268      FGR   1  27612064  27635277
#>   strand        biotype
#> 1     -1 protein_coding
#> 2      1 protein_coding
#> 3     -1 protein_coding
#> 4     -1 protein_coding
#> 5      1 protein_coding
#> 6     -1 protein_coding
#>                                                                                                  description
#> 1                                                          tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
#> 2                                                            tenomodulin [Source:HGNC Symbol;Acc:HGNC:17757]
#> 3 dolichyl-phosphate mannosyltransferase polypeptide 1, catalytic subunit [Source:HGNC Symbol;Acc:HGNC:3005]
#> 4                                               SCY1-like, kinase-like 3 [Source:HGNC Symbol;Acc:HGNC:19285]
#> 5                                    chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
#> 6                          FGR proto-oncogene, Src family tyrosine kinase [Source:HGNC Symbol;Acc:HGNC:3697]
```
