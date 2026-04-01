# Transform Raw Scores to Fractions

Transforms raw xCell scores to estimated cell fractions.

## Usage

``` r
transformScores(scores, fit.vals, scale = TRUE, fn = NULL)
```

## Arguments

- scores:

  Raw scores from rawEnrichmentAnalysis.

- fit.vals:

  Calibration values from spill object.

- scale:

  Logical indicating whether to use scaling. Default is \`TRUE\`.

- fn:

  Character string for saving scores. Default is \`NULL\`.

## Value

Transformed xCell scores matrix.
