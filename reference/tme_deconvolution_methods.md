# TME Deconvolution Methods

A named vector of supported tumor microenvironment (TME) deconvolution
methods in the IOBR package.

## Usage

``` r
tme_deconvolution_methods
```

## Format

A named character vector where names are display names and values are
internal method names.

## Value

Named character vector of available deconvolution methods.

## Details

The methods currently supported are:

- \`mcpcounter\`: MCP-counter for immune and stromal cell populations

- \`epic\`: EPIC for immune, stromal, and cancer cell fractions

- \`xcell\`: xCell for 64 immune and stromal cell types

- \`cibersort\`: CIBERSORT for 22 immune cell types

- \`cibersort_abs\`: CIBERSORT in absolute mode

- \`ips\`: Immunophenoscore calculation

- \`estimate\`: ESTIMATE for stromal/immune/estimate scores

- \`svr\`: Support Vector Regression (custom reference)

- \`lsei\`: Least Squares with Equality/Inequality constraints

- \`timer\`: TIMER for cancer-specific immune estimation

- \`quantiseq\`: quanTIseq for RNA-seq immune deconvolution

## Examples

``` r
# List all available TME deconvolution methods
tme_deconvolution_methods
#>         MCPcounter               EPIC              xCell          CIBERSORT 
#>       "mcpcounter"             "epic"            "xcell"        "cibersort" 
#> CIBERSORT Absolute                IPS           ESTIMATE                SVR 
#>    "cibersort_abs"              "ips"         "estimate"              "svr" 
#>               lsei              TIMER          quanTIseq 
#>             "lsei"            "timer"        "quantiseq" 

# Get internal method name for a specific method
tme_deconvolution_methods["MCPcounter"]
#>   MCPcounter 
#> "mcpcounter" 
```
