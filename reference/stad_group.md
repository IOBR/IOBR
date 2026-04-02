# Example Clinical Data for TCGA-STAD Gastric Cancer Analysis

A small example dataset demonstrating the clinical data structure
required for TCGA Stomach Adenocarcinoma (STAD) analysis using IOBR
package functions. Contains simulated clinical variables, molecular
subtypes, and survival data.

## Usage

``` r
data(stad_group)
```

## Format

A data frame with patient samples as rows and clinical variables as
columns:

- ID:

  Unique patient identifier (TCGA barcode format)

- stage:

  AJCC pathological stage (Stage_II, Stage_III, Stage_IV)

- status:

  Patient vital status (Alive, Dead, NA)

- Lauren:

  Lauren histological classification (Intestinal, Diffuse, Mixed)

- subtype:

  Molecular subtype classification (EBV, GS)

- EBV:

  Epstein-Barr virus status (Positive, Negative)

- time:

  Overall survival time in months

- OS_status:

  Overall survival status (0=alive, 1=dead)

## Examples

``` r
data(stad_group)
head(stad_group)
#>               ID     stage status     Lauren subtype      EBV  time OS_status
#> 33  TCGA-BR-7196  Stage_IV  Alive      Mixed     EBV Positive 22.20         0
#> 89  TCGA-BR-8686 Stage_III   Dead Intestinal     EBV Positive 21.17         1
#> 96  TCGA-BR-A4J4 Stage_III  Alive      Mixed     EBV Positive  0.53         0
#> 18  TCGA-BR-6455  Stage_II   Dead      Mixed     EBV Positive 14.07         1
#> 200 TCGA-FP-7916 Stage_III   Dead    Diffuse     EBV Positive 14.27         1
#> 95  TCGA-BR-A4IV Stage_III   Dead    Diffuse      GS Negitive 28.97         1
```
