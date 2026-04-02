# Add Risk Score to Dataset

Computes a risk score for each observation based on Cox proportional
hazards regression or binary logistic regression. The function fits the
specified model and returns the dataset with an added risk score column.

## Usage

``` r
add_riskscore(
  input,
  family = c("cox", "binary"),
  target = NULL,
  time = NULL,
  status = NULL,
  vars,
  new_var_name = "riskscore"
)
```

## Arguments

- input:

  Data frame containing the variables for analysis.

- family:

  Character string specifying the model family: \`"cox"\` for Cox
  proportional hazards regression or \`"binary"\` for logistic
  regression. Default is \`"cox"\`.

- target:

  Character string specifying the target variable name. Required when
  \`family = "binary"\`.

- time:

  Character string specifying the time-to-event variable name. Required
  when \`family = "cox"\`.

- status:

  Character string specifying the event status variable name. Required
  when \`family = "cox"\`.

- vars:

  Character vector of variable names to include in the model.

- new_var_name:

  Character string specifying the name for the new risk score column.
  Default is \`"riskscore"\`.

## Value

Data frame identical to \`input\` with an additional column containing
risk scores (linear predictors for Cox models or predicted probabilities
for logistic models).

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
input_data <- data.frame(
  time = rexp(100),
  status = rbinom(100, 1, 0.5),
  age = rnorm(100, 60, 10),
  score1 = rnorm(100),
  score2 = rnorm(100)
)
result <- add_riskscore(
  input_data,
  time = "time", status = "status",
  vars = c("age", "score1", "score2")
)
head(result$riskscore)
#> [1]  0.14482695 -0.01998538 -0.16480176  0.02682253 -0.15742721 -0.08938045
```
