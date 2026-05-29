# Get IOBR Cache Directory

Returns the current cache directory for IOBR downloaded data. To comply
with CRAN policies, the default cache directory is a session-specific
temporary directory. Users can opt-in to a persistent cache by setting
\`options(IOBR.cache_dir = "your/path")\` or using
\`set_iobr_cache_dir()\`.

The cache directory is determined in the following priority order: 1.
Function argument \`cache_dir\` (if provided) 2. Option
\`IOBR.cache_dir\` (if set via \`options()\`) 3. Default: A
session-specific temporary directory (\`file.path(tempdir(),
"IOBR_cache")\`)

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
# Get current cache directory (defaults to tempdir)
get_iobr_cache_dir()
#> [1] "/tmp/RtmpVQsgYA"

# Set custom cache directory via options (use tempdir() for examples)
options(IOBR.cache_dir = tempdir())
get_iobr_cache_dir()
#> [1] "/tmp/RtmpVQsgYA"
```
