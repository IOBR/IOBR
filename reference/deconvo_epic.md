# Deconvolve Immune Microenvironment Using EPIC

Estimates immune cell fractions using EPIC algorithm.

## Usage

``` r
deconvo_epic(eset, project = NULL, tumor = TRUE)
```

## Arguments

- eset:

  Gene expression matrix with genes as row names.

- project:

  Optional project name. Default is \`NULL\`.

- tumor:

  Logical indicating tumor (\`TRUE\`) or normal (\`FALSE\`) samples.

## Value

Data frame with EPIC cell fraction estimates. Columns suffixed with
\`\_EPIC\`.

## Author

Dongqiang Zeng

## Examples

``` r
if (FALSE) { # \dontrun{
TRef <- load_data("TRef")
if (!is.null(TRef)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(nrow(TRef$refProfiles) * 2), nrow(TRef$refProfiles), 2)
  rownames(sim_eset) <- rownames(TRef$refProfiles)
  colnames(sim_eset) <- paste0("Sample", 1:2)
  result <- deconvo_epic(eset = sim_eset, project = "Example", tumor = TRUE)
  if (!is.null(result)) head(result)
}
} # }
```
