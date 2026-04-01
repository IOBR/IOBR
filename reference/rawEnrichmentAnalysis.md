# Calculate Raw xCell Enrichment Scores

Returns the raw xCell cell types enrichment scores using ssGSEA.

## Usage

``` r
rawEnrichmentAnalysis(expr, signatures, genes, file.name = NULL)
```

## Arguments

- expr:

  Gene expression data matrix with row names as gene symbols.

- signatures:

  GMT object of signatures.

- genes:

  Character vector of genes to use.

- file.name:

  Character string for saving scores. Default is \`NULL\`.

## Value

Matrix of raw xCell scores.
