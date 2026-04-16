# Construct Mutation Matrices from MAF Data

Builds mutation presence/absence matrices from MAF input (file path or
MAF object). Supports multiple categories: all mutations, SNPs, indels,
and frameshift mutations. When category = "multi", returns a list of
matrices for each category. Compatible with TCGA-formatted data.

## Usage

``` r
make_mut_matrix(
  maf = NULL,
  mut_data = NULL,
  isTCGA = TRUE,
  category = c("multi", "all", "snp", "indel", "frameshift"),
  Tumor_Sample_Barcode = "Tumor_Sample_Barcode",
  Hugo_Symbol = "Hugo_Symbol",
  Variant_Classification = "Variant_Classification",
  Variant_Type = "Variant_Type"
)
```

## Arguments

- maf:

  Character or MAF object. Path to MAF file or an already loaded MAF
  object.

- mut_data:

  Data frame or NULL. Preloaded MAF-like data (used if 'maf' is NULL).

- isTCGA:

  Logical. Whether the MAF follows TCGA conventions. Default is TRUE.

- category:

  Character. Mutation category: "all", "snp", "indel", "frameshift", or
  "multi". Default is "multi".

- Tumor_Sample_Barcode:

  Character. Column name for tumor sample IDs. Default is
  "Tumor_Sample_Barcode".

- Hugo_Symbol:

  Character. Column name for gene symbols. Default is "Hugo_Symbol".

- Variant_Classification:

  Character. Column name for variant classification (e.g.,
  Frame_Shift_Del). Default is "Variant_Classification".

- Variant_Type:

  Character. Column name for variant type (e.g., SNP, INS, DEL). Default
  is "Variant_Type".

## Value

List of mutation matrices (if category = "multi") or a single matrix for
the specified category.

## Note

Some users may encounter errors from upstream data import (e.g. "Can't
combine ..\$Tumor_Seq_Allele2" when using TCGAbiolinks or
TCGAmutations). This is due to inconsistent column types in the source
MAF tables, not an issue of this function. Please ensure your MAF or
merged data frame uses consistent column types (e.g. convert allele
columns to character before input).

## Author

Dongqian Zeng

Shixiang Huang

## Examples

``` r
# \donttest{
# See maftools or TCGAbiolinks documentation for obtaining MAF input
mut_list <- make_mut_matrix(maf = maf, isTCGA = TRUE, category = "multi")
#> Error: object 'maf' not found
# }
```
