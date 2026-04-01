# Split Data into Training and Testing Sets

Divides dataset into training and testing sets using random sampling.
Maintains data integrity for both binomial and survival analysis types.

## Usage

``` r
SplitTrainTest(x, y, train_ratio, type = c("binomial", "survival"), seed)
```

## Arguments

- x:

  Predictor matrix or data frame.

- y:

  Outcome vector (binomial) or matrix with time/status (survival).

- train_ratio:

  Proportion for training (0-1). Default is \`0.7\`.

- type:

  Analysis type: \`"binomial"\` or \`"survival"\`.

- seed:

  Random seed for reproducibility.

## Value

List containing:

- train.x:

  Training predictors matrix

- train.y:

  Training outcomes

- test.x:

  Testing predictors matrix

- test.y:

  Testing outcomes

- train_sample:

  Indices of training samples

## Examples

``` r
# \donttest{
data_matrix <- matrix(rnorm(200), ncol = 2)
outcome_vector <- rbinom(100, 1, 0.5)
split_data <- SplitTrainTest(
  data_matrix, outcome_vector,
  train_ratio = 0.7,
  type = "binomial", seed = 123
)
# }
```
