# Break Time Into Blocks

Divides time duration into specified blocks for analysis.

## Usage

``` r
calculate_break_month(input, block = 6, time_type = c("month", "day"))
```

## Arguments

- input:

  Numeric vector of time durations.

- block:

  Number of blocks. Default is 6.

- time_type:

  Units: "month" or "day". Default is "month".

## Value

Numeric vector of breakpoints, rounded to nearest multiple of 5.

## Author

Dongqiang Zeng

## Examples

``` r
time_data <- c(24, 36, 12, 48)
blocks <- calculate_break_month(input = time_data)
#> ℹ Maximum follow-up time is 48 months; divided into 6 sections
```
