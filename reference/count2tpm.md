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
#> ℹ Found 1679 duplicate symbols. Using "mean" for ranking.
#> ✔ Reduced to 54658 unique genes
head(eset)
#>         TCGA-BR-6455 TCGA-BR-7196 TCGA-BR-8371 TCGA-BR-8380 TCGA-BR-8592
#> MT-CO1     25097.553     24055.97   59199.1348    50164.982     47776.98
#> MT-ND4     14525.308     14963.92   48195.6709    32185.131     36283.39
#> MT-CO3     12465.198     14835.81   47590.7666    34559.375     35808.71
#> IGKC        9897.134     78664.05     713.5372     4372.862     27027.24
#> MT-CO2     12920.617     13063.07   33664.6705    21527.455     25992.41
#> MT-RNR2    18660.657     10908.66   41096.2915    20342.721     24383.75
#>         TCGA-BR-8686 TCGA-BR-A4IV TCGA-BR-A4J4 TCGA-BR-A4J9 TCGA-FP-7916
#> MT-CO1      34875.46     50009.92    27982.230     63380.90     26287.97
#> MT-ND4      31111.09     35885.02    49689.108     39741.24     14473.97
#> MT-CO3      27071.10     37804.39    22623.828     43837.58     17801.27
#> IGKC        50262.87     10617.20     2900.288     21222.09     70074.94
#> MT-CO2      20287.58     32023.19    54797.893     30150.49     13184.96
#> MT-RNR2     25519.06     36721.41    31241.258     31834.37     12436.25
```
