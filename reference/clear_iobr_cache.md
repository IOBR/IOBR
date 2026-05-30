# Clear IOBR Data Cache

Removes all cached data files downloaded from GitHub.

## Usage

``` r
clear_iobr_cache(cache_dir = NULL)
```

## Arguments

- cache_dir:

  Character string. Custom cache directory. If NULL, uses the option
  \`IOBR.cache_dir\` or the default system cache location.

## Value

Invisible NULL. Called for side effects of clearing the cache.

## Examples

``` r
clear_iobr_cache()
#> ✔ Cache cleared: 1 file(s) removed from /tmp/RtmpRGChSj/IOBR_cache
```
