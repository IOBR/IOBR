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
