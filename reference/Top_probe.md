# Top Probe Selector

Extracts the top \`i\` probes based on their ordering in the provided
data frame. If the number of rows is less than or equal to \`i\`,
returns all probes.

## Usage

``` r
Top_probe(dat, i)
```

## Arguments

- dat:

  Data frame containing a column named "probe".

- i:

  Integer. Number of top probes to return.

## Value

Character vector containing the names of the top \`i\` probes.

## Examples

``` r
dat <- data.frame(
  probe = c("Probe1", "Probe2", "Probe3", "Probe4", "Probe5"),
  value = c(5, 3, 2, 4, 1)
)
top_probes <- Top_probe(dat, 3)
print(top_probes)
#> [1] "Probe1" "Probe2" "Probe3"
```
