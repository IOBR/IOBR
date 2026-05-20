# Add Custom Download Mirror

Adds a custom mirror URL to the default mirrors for the current session.
The mirror URL should be a base URL. If it contains 'github.com', it
will be treated as a GitHub proxy and the relative path to IOBR releases
will be appended. Otherwise, it will be treated as a direct repository
path.

## Usage

``` r
add_iobr_mirror(url, position = c("first", "last", "before_github"))
```

## Arguments

- url:

  Character string. The mirror URL to add.

- position:

  Character. Where to add the mirror: "first", "last", or
  "before_github". Default: "first".

## Value

Invisibly returns the updated mirror list.

## Examples

``` r
# \donttest{
# Add a custom mirror to try first
add_iobr_mirror("https://my-mirror.com/https://github.com")
#> ✔ Added mirror <https://my-mirror.com/https://github.com> to position: "first"
#> ℹ Current mirrors: 13 total

# Add mirror to try before default GitHub
add_iobr_mirror("https://fast-mirror.org", position = "before_github")
#> ✔ Added mirror <https://fast-mirror.org> to position: "before_github"
#> ℹ Current mirrors: 14 total

# Download with the new mirror
data <- download_iobr_data("BRef")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "BRef"
# }
```
