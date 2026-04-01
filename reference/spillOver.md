# Adjust Scores Using Spillover Compensation

Adjusts xCell scores using spillover compensation matrix.

## Usage

``` r
spillOver(transformedScores, K, alpha = 0.5, file.name = NULL)
```

## Arguments

- transformedScores:

  Transformed scores from transformScores.

- K:

  Spillover matrix.

- alpha:

  Spillover alpha parameter. Default is \`0.5\`.

- file.name:

  Character string for saving scores. Default is \`NULL\`.

## Value

Adjusted xCell scores matrix.
