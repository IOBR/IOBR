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
if (interactive()) {
  add_iobr_mirror("https://my-mirror.com/https://github.com")
}
```
