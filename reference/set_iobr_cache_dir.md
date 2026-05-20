# Set IOBR Cache Directory

Sets a custom cache directory for IOBR downloaded data. This is useful
when you want to store cached data in a specific location, such as a
shared network drive or a custom directory.

## Usage

``` r
set_iobr_cache_dir(path, create = TRUE)
```

## Arguments

- path:

  Character string. The path to the cache directory.

- create:

  Logical. Whether to create the directory if it doesn't exist. Default:
  TRUE.

## Value

Invisibly returns the cache directory path.

## Examples

``` r
# \donttest{
# Set a custom cache directory (use tempdir() for examples)
set_iobr_cache_dir(tempdir())
#> ✔ IOBR cache directory set to: /tmp/RtmpTveXhQ

# Check the current cache directory
get_iobr_cache_dir()
#> [1] "/tmp/RtmpTveXhQ"

# Download data will now use the custom cache
data <- download_iobr_data("lm22")
#> ℹ Loading cached data: "lm22"
# }
```
