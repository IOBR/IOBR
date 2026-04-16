# Download IOBR Data from GitHub with Mirror Support

Downloads large datasets from GitHub releases to avoid CRAN size limits.
Supports multiple download mirrors for users in different regions. Data
is cached locally after first download.

## Usage

``` r
download_iobr_data(
  name,
  force = FALSE,
  verbose = TRUE,
  mirrors = get_default_mirrors()
)
```

## Arguments

- name:

  Character string. Name of the dataset to download.

- force:

  Logical. Whether to force re-download even if cached. Default: FALSE.

- verbose:

  Logical. Whether to print progress messages. Default: TRUE.

- mirrors:

  Character vector. URLs of mirrors to try. Default uses
  get_default_mirrors().

## Value

The requested dataset.

## Examples

``` r
# \donttest{
# Download TCGA STAD signature data
tcga_sig <- download_iobr_data("tcga_stad_sig")
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "tcga_stad_sig"

# Download with custom mirrors
eset <- download_iobr_data("eset_stad",
  mirrors = c(
    "https://ghproxy.vip/https://github.com",
    "https://gh-proxy.org/https://github.com"
  )
)
#> ℹ Loading cached data: "eset_stad"
# }
```
