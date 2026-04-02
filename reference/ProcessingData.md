# Process Data for Model Construction

Preprocesses data for binomial or survival analysis. Aligns and filters
data based on sample IDs, optionally scales data, and ensures
appropriate data types. Handles missing values by removing columns with
NA values.

## Usage

``` r
ProcessingData(x, y, scale, type = c("binomial", "survival"))
```

## Arguments

- x:

  Data frame of predictors with first column as IDs.

- y:

  Data frame of outcomes with first column as IDs. For survival, expects
  two additional columns for time and status.

- scale:

  Logical indicating whether to scale predictors.

- type:

  Character string: \`"binomial"\` or \`"survival"\`.

## Value

List containing:

- x_scale:

  Processed predictor matrix

- y:

  Processed outcome variable

- x_ID:

  Sample IDs

## Examples

``` r
x <- data.frame(ID = 1:10, predictor1 = rnorm(10), predictor2 = rnorm(10))
y <- data.frame(ID = 1:10, outcome = sample(c(0, 1), 10, replace = TRUE))
result <- ProcessingData(x, y, scale = TRUE, type = "binomial")
#> ℹ Converting outcome to factor
```
