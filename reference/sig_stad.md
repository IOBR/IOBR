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
