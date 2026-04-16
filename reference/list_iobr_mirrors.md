# List Current Download Mirrors

Returns the current list of download mirrors.

## Usage

``` r
list_iobr_mirrors()
```

## Value

Character vector of mirror URLs.

## Examples

``` r
list_iobr_mirrors()
#> Current IOBR download mirrors:
#> ==============================
#> 1. https://my-mirror.com/https://github.com
#> 2. https://fast-mirror.org
#> 3. https://github.com
#> 4. https://ghproxy.vip/https://github.com
#> 5. https://gh-proxy.org/https://github.com
#> 6. https://ghfast.top/https://github.com
#> ==============================
#> Total: 6 mirrors
```
