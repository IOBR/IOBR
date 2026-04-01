# Get Outlier Genes

Identifies outlier genes from multiple cancer datasets. Treats the top 5
expressed genes in each sample as outliers and returns unique outlier
genes.

## Usage

``` r
GetOutlierGenes(cancers)
```

## Arguments

- cancers:

  Data frame. One column containing paths to gene expression files.

## Value

Vector of unique gene names identified as outliers.

## Examples

``` r
tf1 <- tempfile(fileext = ".tsv")
tf2 <- tempfile(fileext = ".tsv")

expr1 <- data.frame(
  gene = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"),
  Sample1 = c(10, 50, 30, 80, 60, 20),
  Sample2 = c(15, 40, 25, 90, 55, 10)
)

expr2 <- data.frame(
  gene = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"),
  Sample3 = c(100, 20, 10, 60, 30, 80),
  Sample4 = c(95, 25, 15, 70, 35, 85)
)

write.table(expr1, tf1, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(expr2, tf2, sep = "\t", row.names = FALSE, quote = FALSE)

cancers <- data.frame(ExpressionFiles = c(tf1, tf2))
outlier_genes <- GetOutlierGenes(cancers)
```
