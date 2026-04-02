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

## Examples

``` r
data(eset_stad)
head(eset_stad)
#>                 TCGA-BR-6455 TCGA-BR-7196 TCGA-BR-8371 TCGA-BR-8380
#> ENSG00000000003         8006         2114          767         1556
#> ENSG00000000005            1            0            5            5
#> ENSG00000000419         3831         2600         1729         1760
#> ENSG00000000457         1126          745         1040         1260
#> ENSG00000000460          857          463          231          432
#> ENSG00000000938          758         1126          557          557
#>                 TCGA-BR-8592 TCGA-BR-8686 TCGA-BR-A4IV TCGA-BR-A4J4
#> ENSG00000000003         2806         2923         1524         7208
#> ENSG00000000005           60            1           22            2
#> ENSG00000000419         2273         1934         2838         4418
#> ENSG00000000457         1814          707         1683         1335
#> ENSG00000000460          635          323          270          423
#> ENSG00000000938          828          666          760          597
#>                 TCGA-BR-A4J9 TCGA-FP-7916
#> ENSG00000000003          711         2747
#> ENSG00000000005            0            3
#> ENSG00000000419         2426         2824
#> ENSG00000000457         1590         1672
#> ENSG00000000460          276          773
#> ENSG00000000938          370          688
```
