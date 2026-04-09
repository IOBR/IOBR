# Add Custom Download Mirror

Adds a custom mirror URL to the default mirrors for the current session.
The mirror URL should be a base URL that will be prepended to GitHub
paths.

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
if (FALSE) { # \dontrun{
# Add a custom mirror to try first
add_iobr_mirror("https://my-mirror.com/https://github.com")

# Add mirror to try before default GitHub
add_iobr_mirror("https://fast-mirror.org", position = "before_github")

# Download with the new mirror
data <- download_iobr_data("BRef")
} # }
```
