# Affymetrix Human Genome U133 Plus 2.0 Array Annotation

Probe annotation for the Affymetrix Human Genome U133 Plus 2.0
microarray. Provides mapping between probe IDs and gene symbols for
microarray data analysis.

## Usage

``` r
data(anno_hug133plus2)
```

## Format

A data frame with columns:

- probe_id:

  Affymetrix probe identifier.

- symbol:

  Official gene symbol.

## Source

Affymetrix annotation files (HuGene-1_0-st-v1)

## Examples

``` r
data(anno_hug133plus2)
head(anno_hug133plus2)
#> # A tibble: 6 × 2
#>   probe_id  symbol 
#>   <fct>     <fct>  
#> 1 1007_s_at MIR4640
#> 2 1053_at   RFC2   
#> 3 117_at    HSPA6  
#> 4 121_at    PAX8   
#> 5 1255_g_at GUCA1A 
#> 6 1294_at   MIR5193
```
