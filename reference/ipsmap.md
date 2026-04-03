# Map Score to Immunophenoscore

Maps input score to Immunophenoscore (IPS) on a 0-10 scale. Scores
\\\le\\0 map to 0, scores \\\ge\\3 map to 10, and intermediate scores
are linearly scaled.

## Usage

``` r
ipsmap(x)
```

## Arguments

- x:

  Numeric value representing the aggregate z-score.

## Value

Integer value between 0 and 10 representing the Immunophenoscore.

## Examples

``` r
ips <- ipsmap(2.5)
ips <- ipsmap(-1)
ips <- ipsmap(5)
```
