# Illumina Microarray Annotation

Probe annotation for Illumina microarray platforms. Provides mapping
between Illumina probe IDs and gene symbols for microarray data
analysis.

## Usage

``` r
data(anno_illumina)
```

## Format

A data frame with columns:

- probe_id:

  Illumina probe identifier.

- symbol:

  Official gene symbol.

## Source

Illumina annotation manifest files

## Examples

``` r
data(anno_illumina)
head(anno_illumina)
#> # A tibble: 6 × 2
#>   probe_id     symbol   
#>   <fct>        <fct>    
#> 1 ILMN_1725881 LOC23117 
#> 2 ILMN_1910180 HS.575038
#> 3 ILMN_1804174 FCGR2B   
#> 4 ILMN_1796063 TRIM44   
#> 5 ILMN_1811966 LOC653895
#> 6 ILMN_1668162 DGAT2L3  
```
