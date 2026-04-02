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

## Examples

``` r
data(tcga_stad_pdata)
head(tcga_stad_pdata)
#> # A tibble: 6 × 17
#>   ID        stage status Lauren subtype EBV   TMEscore_plus TMEscore_plus_binary
#>   <fct>     <fct> <fct>  <fct>  <fct>   <chr>         <dbl> <fct>               
#> 1 TCGA-3M-… Stag… Alive  Mixed  NA      NE           0.0874 Low                 
#> 2 TCGA-B7-… Stag… Alive  Diffu… EBV     Posi…        2.03   High                
#> 3 TCGA-B7-… Stag… Alive  Diffu… NA      NE          -0.0423 Low                 
#> 4 TCGA-B7-… Stag… Alive  Intes… NA      NE          -0.686  Low                 
#> 5 TCGA-B7-… Stag… Alive  Intes… NA      NE           2.51   High                
#> 6 TCGA-B7-… Stag… Alive  Mixed  NA      NE          -0.627  Low                 
#> # ℹ 9 more variables: time <dbl>, OS_status <int>, ARID1A <chr>, PIK3CA <chr>,
#> #   MALAT1 <dbl>, GZMB <dbl>, GNLY <dbl>, CD274 <dbl>, HOTAIR <dbl>
```
