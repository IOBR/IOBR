# Extract Best Cutoff and Add Binary Variable to Data Frame

Determines the optimal cutoff point for a continuous variable in
survival analysis using the maximally selected rank statistics method.
Creates a binary variable based on the identified cutoff and adds it to
the input data frame.

## Usage

``` r
best_cutoff(
  pdata,
  variable,
  time = "time",
  status = "status",
  print_result = TRUE
)
```

## Arguments

- pdata:

  Data frame containing survival information and the continuous
  variable.

- variable:

  Character string specifying the name of the continuous variable for
  which the optimal cutoff should be determined.

- time:

  Character string specifying the column name containing time-to-event
  data. Default is \`"time"\`.

- status:

  Character string specifying the column name containing event status
  (censoring information). Default is \`"status"\`.

- print_result:

  Logical indicating whether to print detailed results including cutoff
  value and Cox model summaries. Default is \`TRUE\`.

## Value

Data frame identical to \`pdata\` with an additional binary column named
\`\<variable\>\_binary\` containing "High" and "Low" categories based on
the optimal cutoff.

## Author

Dongqiang Zeng

## Examples

``` r
set.seed(123)
pdata <- data.frame(
  time = rexp(100),
  status = rbinom(100, 1, 0.5),
  score = rnorm(100, mean = 50, sd = 10)
)
result <- best_cutoff(pdata, variable = "score", print_result = FALSE)
#> ✔ Best cutoff for "score": 46.503
table(result$score_binary)
#> 
#>  Low High 
#>   45   55 
```
