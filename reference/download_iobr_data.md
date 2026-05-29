# Download IOBR Data from GitHub with Mirror Support

Downloads large datasets from GitHub releases to avoid CRAN size limits.
Supports multiple download mirrors for users in different regions. Data
is cached locally after first download. Cache directory can be
customized via \`options(IOBR.cache_dir = "your/path")\`.

## Usage

``` r
download_iobr_data(
  name,
  force = FALSE,
  verbose = TRUE,
  mirrors = get_default_mirrors(),
  cache_dir = NULL
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

- cache_dir:

  Character string. Custom cache directory. If NULL, uses the option
  \`IOBR.cache_dir\` or the default system cache location.

## Value

The requested dataset.

## Examples

``` r
if (interactive()) {
  tcga_sig <- download_iobr_data("tcga_stad_sig")
}
```
