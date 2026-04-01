# Toy STAD Phenotype Data

A data frame containing clinical and pathological annotations for the
TCGA stomach adenocarcinoma (STAD) cohort. Each row corresponds to one
tumour sample and can be matched to the columns of `eset_stad` via the
`ID` column. This dataset is typically used together with `eset_stad` in
examples of survival analysis, subgroup comparison and immune
deconvolution in the IOBR package.

## Usage

``` r
data(pdata_stad)
```

## Format

A data frame with one row per TCGA-STAD sample and 8 variables:

- ID:

  Character. TCGA sample barcode, matching the column names of
  `eset_stad`.

- stage:

  Factor. Pathological stage (e.g. `"Stage_I"`, `"Stage_II"`,
  `"Stage_III"`, `"Stage_IV"`).

- status:

  Factor. Vital status at last follow-up (`"Alive"` or `"Dead"`).

- Lauren:

  Factor. Lauren classification of gastric cancer (`"Intestinal"`,
  `"Diffuse"`, `"Mixed"` or `NA`).

- subtype:

  Factor. Molecular subtype (e.g. `"CIN"`, `"EBV"`, `"GS"`, `"MSI"`).

- EBV:

  Factor. EBV status of the tumour (`"Positive"` or `"Negative"`).

- time:

  Numeric. Overall survival or follow-up time, typically measured in
  months.

- OS_status:

  Integer/binary. Overall survival status indicator. (`1` = death, `0` =
  censored)
