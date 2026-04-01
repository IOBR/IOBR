# TCGA-STAD Clinical and Molecular Annotation Data

Clinical, molecular, and signature score data for TCGA stomach
adenocarcinoma (STAD) samples. Includes patient demographics, tumor
characteristics, molecular subtypes, and pre-computed signature scores.

## Usage

``` r
data(tcga_stad_pdata)
```

## Format

A data frame with samples as rows and variables as columns:

- ID – TCGA sample barcode

- stage – tumor stage (Stage_I, Stage_II, etc.)

- status – vital status (Alive/Dead)

- Lauren – histological classification

- subtype – molecular subtype (EBV, MSI, GS, CN, ...)

- EBV – EBV infection status

- TMEscore_plus – continuous tumor-micro-environment score

- TMEscore_plus_binary – High/Low TME classification

- time – follow-up time (months)

- OS_status – 0 = alive, 1 = dead

- ARID1A, PIK3CA – driver-gene mutation status

- MALAT1 – lncRNA expression

- remaining columns – gene-expression values and additional
  clinical/molecular annotations

## Source

The Cancer Genome Atlas Stomach Adenocarcinoma (TCGA-STAD)

## References

Cancer Genome Atlas Research Network. Comprehensive molecular
characterization of gastric adenocarcinoma. Nature 513, 202-209 (2014).
doi:10.1038/nature13480
