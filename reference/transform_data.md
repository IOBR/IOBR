# Transform NA, Inf, or Zero Values in Data

Replaces NA, Inf, or zero values in specified columns of a data frame
with a user-defined value or the column mean.

## Usage

``` r
transform_data(data, feature, data_type = c("NA", "Inf", "zero"), into = 0)
```

## Arguments

- data:

  Data frame. Input data to be transformed.

- feature:

  Character vector. Column names in \`data\` to apply transformation.

- data_type:

  Character. Type of value to replace: \`"NA"\`, \`"Inf"\`, or
  \`"zero"\`.

- into:

  Value to replace specified type with. Default is 0. If \`"mean"\`,
  replaces with column mean (excluding NA/Inf values as appropriate).

## Value

Data frame with specified transformations applied to selected features.

## Author

Dongqiang Zeng

## Examples

``` r
data_matrix <- data.frame(
  A = c(1, 2, NA, 4, Inf),
  B = c(Inf, 2, 3, 4, 5),
  C = c(0, 0, 0, 1, 2)
)

# Replace NAs with 0
transform_data(data_matrix, feature = c("A", "B"), data_type = "NA")
#>     A   B C
#> 1   1 Inf 0
#> 2   2   2 0
#> 3   0   3 0
#> 4   4   4 1
#> 5 Inf   5 2

# Replace Inf values with the mean of the column
transform_data(data_matrix,
  feature = c("A", "B"),
  data_type = "Inf", into = "mean"
)
#>          A   B C
#> 1 1.000000 3.5 0
#> 2 2.000000 2.0 0
#> 3       NA 3.0 0
#> 4 4.000000 4.0 1
#> 5 2.333333 5.0 2

# Replace zeros with -1 in column C
transform_data(data_matrix, feature = "C", data_type = "zero", into = -1)
#>     A   B  C
#> 1   1 Inf -1
#> 2   2   2 -1
#> 3  NA   3 -1
#> 4   4   4  1
#> 5 Inf   5  2
```
