# Annotate Gene Expression Matrix and Remove Duplicated Genes

Annotates an expression matrix with gene symbols using provided
annotation data, filters out missing or invalid symbols, handles
duplicate gene entries, and removes uninformative rows. The function
supports multiple aggregation methods for resolving duplicate gene
symbols.

## Usage

``` r
anno_eset(
  eset,
  annotation,
  symbol = "symbol",
  probe = "probe_id",
  method = "mean"
)
```

## Arguments

- eset:

  Expression matrix or ExpressionSet object containing gene expression
  data.

- annotation:

  Data frame containing annotation information for probes. Built-in
  options include \`anno_hug133plus2\`, \`anno_rnaseq\`, and
  \`anno_illumina\`.

- symbol:

  Character string specifying the column name in \`annotation\` that
  represents gene symbols. Default is \`"symbol"\`.

- probe:

  Character string specifying the column name in \`annotation\` that
  represents probe identifiers. Default is \`"probe_id"\`.

- method:

  Character string specifying the aggregation method for duplicate gene
  symbols. Options are \`"mean"\`, \`"sum"\`, or \`"sd"\`. Default is
  \`"mean"\`.

## Value

Annotated and cleaned gene expression matrix with gene symbols as row
names.

## Details

The function performs the following operations:

1.  Filters probes with missing symbols or labeled as \`"NA_NA"\`

2.  Matches probes between expression set and annotation data

3.  Merges annotation with expression data

4.  Handles duplicate gene symbols using specified aggregation method

5.  Removes rows with all zeros, all NAs, or missing values in the first
    column

## Author

Dongqiang Zeng

## Examples

``` r
# Annotate Affymetrix microarray data
eset_gse62254 <- load_data("eset_gse62254")
anno_hug133plus2 <- load_data("anno_hug133plus2")
eset <- anno_eset(eset = eset_gse62254, annotation = anno_hug133plus2)
#> ℹ Row number of original eset: 54675
#> ✔ 83% of probes in expression set were annotated
#> ℹ Found 23366 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 21752
head(eset)
#>              GSM1523727 GSM1523728 GSM1523729
#> SH3KBP1        4.327974   4.316195   4.351425
#> EEF1A1         4.293762   4.291038   4.262199
#> COX2           4.250288   4.283714   4.270508
#> RPL41          4.246149   4.246808   4.257940
#> ND4            4.285322   4.218556   4.244270
#> LOC101928826   4.219303   4.219670   4.213252

# Annotate RNA-seq data with Ensembl IDs
eset_stad <- load_data("eset_stad")
anno_grch38 <- load_data("anno_grch38")
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
head(eset)
#>         TCGA-BR-6455 TCGA-BR-7196 TCGA-BR-8371 TCGA-BR-8380 TCGA-BR-8592
#> MT-CO1        849866       858644      2296926      1851900      2317798
#> MT-ND4        439551       477310      1671108      1061786      1573003
#> MYH11          49090       231644      1891538      2278901      2268508
#> MT-RNR2       638863       393662      1612115       759255      1195967
#> FLNA           78921       291182      1255093      1469937      1586939
#> ACTB          367641       435945       806904       864952       872777
#>         TCGA-BR-8686 TCGA-BR-A4IV TCGA-BR-A4J4 TCGA-BR-A4J9 TCGA-FP-7916
#> MT-CO1        713603      2547162      1172210      3164451      1019632
#> MT-ND4        568875      1633347      1860155      1773153       501694
#> MYH11          89832      1867340        58720      2473820       103143
#> MT-RNR2       527914      1890956      1323163      1606934       487683
#> FLNA          105371      1483645       129692      1622593       190714
#> ACTB          200980      1231785       563329      1041598       540407
```
