# GSVA API Version Detection

Detects whether the installed GSVA package supports the new
parameter-based API (gsvaParam/ssgseaParam) or the old direct argument
API. This function is used internally to ensure compatibility across
different GSVA package versions (1.50.0+ vs older versions).

## Usage

``` r
gsva_use_new_api()
```

## Value

A list with two elements: - \`use_new_api\`: Logical indicating whether
to use the new API (\`TRUE\`) or old API (\`FALSE\`) - \`gsva_version\`:
Character string of the installed GSVA version, or \`"not installed"\`
if not available

## Examples

``` r
# Detect GSVA API version (only runs if GSVA is installed)
api_info <- IOBR:::gsva_use_new_api()
print(api_info)
#> $use_new_api
#> [1] TRUE
#> 
#> $gsva_version
#> [1] "2.4.8"
#> 
```
