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
#> 4. https://gh-proxy.com/https://github.com
#> 5. https://ghproxy.net/https://github.com
#> 6. https://moeyy.cn/gh-proxy/https://github.com
#> 7. https://github.akams.cn/https://github.com
#> 8. http://toolwa.com/github/https://github.com
#> 9. https://v6.gh-proxy.org/https://github.com
#> 10. https://gh-proxy.org/https://github.com
#> 11. https://ghfast.top/https://github.com
#> 12. https://download.githubcdn.com?url=https://github.com
#> 13. https://proxy.gitwarp.top/https://github.com
#> ==============================
#> Total: 13 mirrors
```
