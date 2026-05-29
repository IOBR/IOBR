# Convert Read Counts to Transcripts Per Million (TPM)

Transforms gene expression count data into Transcripts Per Million (TPM)
values, normalizing for gene length and library size. Supports multiple
gene ID types and can retrieve gene length information from BioMart or
use local datasets.

## Usage

``` r
count2tpm(
  countMat,
  idType = "Ensembl",
  org = c("hsa", "mmus"),
  source = c("local", "biomart"),
  effLength = NULL,
  id = "id",
  gene_symbol = "symbol",
  length = "eff_length",
  check_data = FALSE
)
```

## Arguments

- countMat:

  Numeric matrix of raw read counts with genes in rows and samples in
  columns.

- idType:

  Character string specifying the gene identifier type. Options are
  \`"Ensembl"\`, \`"Entrez"\`, or \`"Symbol"\`. Default is
  \`"Ensembl"\`.

- org:

  Character string specifying the organism. Options include \`"hsa"\`
  (human) or \`"mmus"\` (mouse). Default is \`"hsa"\`.

- source:

  Character string specifying the source for gene length information.
  Options are \`"biomart"\` (retrieve from Ensembl BioMart) or
  \`"local"\` (use local dataset). Default is \`"local"\`.

- effLength:

  Data frame containing effective gene length information. If \`NULL\`,
  lengths are retrieved based on \`source\`. Default is \`NULL\`.

- id:

  Character string specifying the column name in \`effLength\`
  containing gene identifiers. Default is \`"id"\`.

- gene_symbol:

  Character string specifying the column name in \`effLength\`
  containing gene symbols. Default is \`"symbol"\`.

- length:

  Character string specifying the column name in \`effLength\`
  containing gene lengths. Default is \`"eff_length"\`.

- check_data:

  Logical indicating whether to check for missing values in the count
  matrix. Default is \`FALSE\`.

## Value

Data frame of TPM-normalized expression values with genes in rows and
samples in columns. Gene identifiers are converted to gene symbols in
the output, regardless of the input \`idType\`.

## Author

Wubing Zhang, Dongqiang Zeng, Yiran Fang

## Examples

``` r
# Simulated count data
set.seed(123)
count_matrix <- matrix(
  rpois(100, lambda = 10),
  ncol = 5,
  dimnames = list(
    paste0("ENSG00000", 101:120),
    paste0("Sample", 1:5)
  )
)
# Simulated effective length data
eff_len <- data.frame(
  id = paste0("ENSG00000", 101:120),
  symbol = paste0("Gene", 1:20),
  eff_length = runif(20, 500, 5000)
)
# Transform to TPM using local gene annotation
eset <- count2tpm(
  countMat = count_matrix, effLength = eff_len,
  id = "id", length = "eff_length", gene_symbol = "symbol"
)
#> ℹ No duplicate gene symbols found.
head(eset)
#>         Sample1   Sample2   Sample3   Sample4   Sample5
#> Gene1  17872.27  15587.34  18209.61  24387.20  27170.38
#> Gene2  25679.75  26544.19  23257.31  33978.89  34701.99
#> Gene3  42051.04  13971.40  36724.09  26826.92  30910.36
#> Gene4 146098.96 101936.44 133970.97 130487.57 136681.51
#> Gene5  47308.28  33008.03  48201.25  32863.57  35407.06
#> Gene6  61007.04  37836.40  37295.14  28253.11  53269.60
```
