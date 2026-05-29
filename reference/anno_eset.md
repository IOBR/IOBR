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
# Create a small example expression matrix
eset_mat <- matrix(runif(100), nrow = 10, ncol = 10)
rownames(eset_mat) <- paste0("Probe", 1:10)
colnames(eset_mat) <- paste0("Sample", 1:10)

# Create a matching annotation data frame
anno_df <- data.frame(
  probe_id = paste0("Probe", 1:10),
  symbol = c("Gene1", "Gene1", "Gene2", "Gene3", "Gene4",
             "Gene5", "Gene6", "Gene7", "Gene8", "Gene9")
)

# Annotate
result <- anno_eset(eset = eset_mat, annotation = anno_df)
#> ℹ Row number of original eset: 10
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 1 duplicate symbol, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 9
head(result)
#>         Sample1   Sample2    Sample3     Sample4   Sample5   Sample6   Sample7
#> Gene1 0.7619853 0.6753624 0.66960355 0.852364600 0.8665602 0.3637427 0.7706230
#> Gene2 0.8609894 0.7453418 0.81146366 0.668323644 0.3138426 0.6979466 0.9576697
#> Gene3 0.6735550 0.9425989 0.75916680 0.511314597 0.9594658 0.6844865 0.1586884
#> Gene8 0.1072947 0.5656543 0.92369921 0.003896343 0.8083863 0.8868626 0.9758897
#> Gene6 0.8917137 0.2594267 0.05088011 0.820474512 0.3839367 0.1372436 0.8697061
#> Gene5 0.6931989 0.2977231 0.38090370 0.903362288 0.5314093 0.5546818 0.8731514
#>          Sample8    Sample9  Sample10
#> Gene1 0.38917032 0.95657976 0.6232614
#> Gene2 0.09292582 0.37215774 0.9079464
#> Gene3 0.16180921 0.96261534 0.5694433
#> Gene8 0.56028054 0.50815808 0.7819689
#> Gene6 0.41525745 0.40994593 0.7620283
#> Gene5 0.34181444 0.06126261 0.1168276
```
