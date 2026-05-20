# Get IOBR Cache Directory

Returns the current cache directory for IOBR downloaded data. The cache
directory is determined in the following priority order: 1. Function
argument \`cache_dir\` (if provided) 2. Option \`IOBR.cache_dir\` (if
set via \`options()\`) 3. Default system cache location via
\`tools::R_user_dir()\`

## Usage

``` r
get_iobr_cache_dir(cache_dir = NULL)
```

## Arguments

- cache_dir:

  Optional character string to override the current setting.

## Value

Character string with the cache directory path.

## Examples

``` r
# Get current cache directory
get_iobr_cache_dir()
#> [1] "/tmp/RtmpTveXhQ"

# Set custom cache directory via options (use tempdir() for examples)
options(IOBR.cache_dir = tempdir())
get_iobr_cache_dir()
#> [1] "/tmp/RtmpTveXhQ"
```
