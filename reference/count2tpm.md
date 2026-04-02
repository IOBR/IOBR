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
samples in columns.

## Author

Wubing Zhang, Dongqiang Zeng, Yiran Fang

## Examples

``` r
# Load TCGA count data
eset_stad <- load_data("eset_stad")

# Transform to TPM using local gene annotation
eset <- count2tpm(countMat = eset_stad, source = "local", idType = "ensembl")
#> ℹ Using local annotation (anno_grch38) for TPM conversion
#> ! Omitting 3985 genes without length information
#> Warning: longer object length is not a multiple of shorter object length
#> ℹ No duplicate gene symbols found.
head(eset)
#>                 TCGA-BR-6455 TCGA-BR-7196 TCGA-BR-8371 TCGA-BR-8380
#> ENSG00000000003  46.05359551    5.8697325  1.701672411  6.131800804
#> ENSG00000000005   0.01767806    0.0000000  0.005135287  0.004749853
#> ENSG00000000419  10.77636802   11.5572624  2.789726562  1.253107664
#> ENSG00000000457   4.26856229    0.9786535  3.175097901  3.116723054
#> ENSG00000000460   3.74565258    0.5428787  0.202426512  0.553130705
#> ENSG00000000938   5.84812443    5.9272704  0.413138183  0.517957128
#>                 TCGA-BR-8592 TCGA-BR-8686 TCGA-BR-A4IV TCGA-BR-A4J4
#> ENSG00000000003    8.4607297  8.288326234   0.70255089 16.294383776
#> ENSG00000000005    0.1278757  0.002879186   0.01314326  0.006518611
#> ENSG00000000419    0.9238454  5.373327546   6.26164134  3.198678776
#> ENSG00000000457    0.5530725  1.411757238   4.56153291  1.305919230
#> ENSG00000000460    0.3645667  1.616375779   0.40655373  0.930441480
#> ENSG00000000938    0.7235548  0.663481389   0.51788753  0.472444285
#>                 TCGA-BR-A4J9 TCGA-FP-7916
#> ENSG00000000003    0.3131897   3.94017720
#> ENSG00000000005    0.0000000   0.00251337
#> ENSG00000000419    1.3717145   1.28042954
#> ENSG00000000457    2.0238365   2.94567403
#> ENSG00000000460    0.2368734   5.09922078
#> ENSG00000000938    0.5413752   1.18625545
```
