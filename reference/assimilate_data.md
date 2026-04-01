# Harmonize Two Data Frames by Column Structure

Adds missing columns (filled with \`NA\`) to a secondary data frame so
that its column set and order match a reference data frame. This is
useful when combining data frames from different sources that should
have the same structure but may be missing some columns.

## Usage

``` r
assimilate_data(data_a, data_b)
```

## Arguments

- data_a:

  Data frame. Reference data frame whose column structure should be
  matched.

- data_b:

  Data frame. Data frame to be conformed to \`data_a\`.

## Value

Data frame \`data_b\` with added missing columns (NA-filled) and
reordered to match \`data_a\`.

## Examples

``` r
# Create reference data frame
pdata_a <- data.frame(
  A = 1:5, B = 2:6, C = 3:7, D = 4:8, E = 5:9
)

# Create data frame with subset of columns
pdata_b <- data.frame(A = 1:3, C = 4:6, E = 7:9)

# Harmonize pdata_b to match pdata_a structure
pdata_b_harmonized <- assimilate_data(data_a = pdata_a, data_b = pdata_b)
#> ℹ Adding 2 missing columns: "B" and "D"
print(names(pdata_b_harmonized)) # Now has A, B, C, D, E
#> [1] "A" "B" "C" "D" "E"
```
