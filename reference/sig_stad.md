# TCGA-STAD Gastric Cancer Cohort with Molecular and Clinical Data

A comprehensive dataset containing clinical, pathological, and molecular
data from The Cancer Genome Atlas (TCGA) Stomach Adenocarcinoma (STAD)
project. Includes RNA-seq expression data, survival outcomes, and
pathological features for gastric cancer patients.

## Usage

``` r
data(sig_stad)
```

## Format

A data frame with 374 rows (patients) and 323 variables:

- ID – TCGA barcode

- ProjectID – "TCGA-STAD"

- Technology / platform – sequencing details

- Gender – M/F

- Age – age at diagnosis (years)

- Survival outcomes – RFS_status, OS_time, OS_status

- Pathology – Lauren type, differentiation, AJCC_stage, T/N/M_stage

- 308 omics columns – gene-level RNA-seq counts / TPM / signatures

## Source

The Cancer Genome Atlas (TCGA) Research Network
<https://www.cancer.gov/tcga>

## References

Cancer Genome Atlas Research Network. Comprehensive molecular
characterization of gastric adenocarcinoma. Nature 513, 202-209 (2014).
doi:10.1038/nature13480

Liu J et al. Integrated omics analysis of gastric cancer. Cell Reports
29, 1-15 (2019). doi:10.1016/j.celrep.2019.09.045

## Examples

``` r
data(sig_stad)
head(sig_stad)
#> # A tibble: 6 × 323
#>   ID      ProjectID Technology platform Gender   Age RFS_time RFS_status OS_time
#>   <fct>   <fct>     <fct>      <fct>    <fct>  <dbl>    <dbl>      <int>   <dbl>
#> 1 TCGA-3… TCGA-STAD RNAseq     Illumin… M         70       NA         NA    58.8
#> 2 TCGA-3… TCGA-STAD RNAseq     Illumin… M         51       NA         NA    NA  
#> 3 TCGA-B… TCGA-STAD RNAseq     Illumin… M         62       NA         NA    11.9
#> 4 TCGA-B… TCGA-STAD RNAseq     Illumin… M         52       NA         NA    19.8
#> 5 TCGA-B… TCGA-STAD RNAseq     Illumin… M         74       NA         NA    11.2
#> 6 TCGA-B… TCGA-STAD RNAseq     Illumin… M         51       NA         NA     9.6
#> # ℹ 314 more variables: OS_status <int>, Lauren <fct>, Differtiation <fct>,
#> #   AJCC_Stage_confuse <int>, T_stage <int>, N_stage <int>, M_stage <int>,
#> #   Lymph_node_examined <int>, Positive_lymph_nodes <int>,
#> #   Revisedlocation <fct>, MSI <int>, EBV <int>, Hpylori <int>, Subtype <fct>,
#> #   TP53mutated <int>, B.cells.naive <dbl>, B.cells.memory <dbl>,
#> #   Plasma.cells <dbl>, T.cells.CD8 <dbl>, T.cells.CD4.naive <dbl>,
#> #   T.cells.CD4.memory.resting <dbl>, T.cells.CD4.memory.activated <dbl>, …
```
