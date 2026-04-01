# Format Input Signatures from MSigDB

Reads a GMT file from MSigDB and converts it into a named list of gene
sets suitable for IOBR functions.

## Usage

``` r
format_msigdb(gmt, ont = "term", gene = "gene")
```

## Arguments

- gmt:

  Character string. Path to a GMT file.

- ont:

  Character string. Name of the signature/set column in the parsed GMT
  table. Default is \`"term"\`.

- gene:

  Character string. Name of the gene column in the parsed GMT table.
  Default is \`"gene"\`.

## Value

Named list of character vectors, where each element contains the genes
belonging to one signature.

## Examples

``` r
# \donttest{
tf <- tempfile(fileext = ".gmt")
writeLines(
  c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB\tNA\tTNF\tNFKB1\tNFKB2",
    "HALLMARK_P53_PATHWAY\tNA\tTP53\tMDM2\tCDKN1A"
  ),
  con = tf
)

sig_list <- format_msigdb(tf, ont = "term", gene = "gene")
names(sig_list)
#> [1] "HALLMARK_TNFA_SIGNALING_VIA_NFKB" "HALLMARK_P53_PATHWAY"            
sig_list[[1]]
#> [1] "TNF"   "NFKB1" "NFKB2"
# }
```
