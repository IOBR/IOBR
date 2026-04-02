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

## Examples

``` r
data(eset_blca)
head(eset_blca)
#>                 TCGA-2F-A9KO TCGA-2F-A9KP TCGA-2F-A9KQ TCGA-2F-A9KR
#> ENSG00000000003         6092        11652         5426         4383
#> ENSG00000000005            0            4            1            1
#> ENSG00000000419         3072         2656         1983         2061
#> ENSG00000000457         1302          984         1134         1092
#> ENSG00000000460          779          924          421          386
#> ENSG00000000938          436          116          312          590
#>                 TCGA-2F-A9KT
#> ENSG00000000003         3334
#> ENSG00000000005            0
#> ENSG00000000419         2930
#> ENSG00000000457          496
#> ENSG00000000460          318
#> ENSG00000000938          362
```
