# TIMER Available Cancer Types

Character vector of cancer types supported by TIMER deconvolution. TIMER
signatures are cancer-specific.

## Usage

``` r
timer_available_cancers
```

## Value

Character vector of available cancer type abbreviations.

## Examples

``` r
# List all available cancer types for TIMER
timer_available_cancers
#>  [1] "kich" "blca" "brca" "cesc" "gbm"  "hnsc" "kirp" "lgg"  "lihc" "luad"
#> [11] "lusc" "prad" "sarc" "pcpg" "paad" "tgct" "ucec" "ov"   "skcm" "dlbc"
#> [21] "kirc" "acc"  "meso" "thca" "uvm"  "ucs"  "thym" "esca" "stad" "read"
#> [31] "coad" "chol"

# Check if a cancer type is supported
"brca" %in% timer_available_cancers
#> [1] TRUE
```
