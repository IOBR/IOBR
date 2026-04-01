# Default Pattern List for Name Cleaning

A character vector of common substrings to remove from feature names.
Used in \[remove_names()\] and other helper functions.

## Usage

``` r
patterns_to_na
```

## Format

A character vector of length 12.

## Value

Character vector of patterns to remove.

## Examples

``` r
# View default patterns
patterns_to_na
#>  [1] "_CIBERSORT"  "xCell"       "_EPIC"       "_TIMER"      "_quantiseq" 
#>  [6] "_MCPcounter" "HALLMARK_"   "_cibersort"  "_xcell"      "_epic"      
#> [11] "_timer"      "quanTIseq"   "_mcpcounter" "_hallmark"   "hallmark1"  
#> [16] "hallmark2"   "hallmark3"  
```
