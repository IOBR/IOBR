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
# \donttest{
sig_stad <- load_data("sig_stad")
sig_stad2 <- best_cutoff(
  pdata = sig_stad,
  variable = "TMEscore_CIR",
  time = "OS_time",
  status = "OS_status",
  print_result = TRUE
)
#> ✔ Best cutoff for "TMEscore_CIR": 5.697
#> $best_cutoff
#> [1] "The best cutoff is = 5.697"
#> 
#> $cox_continuous_object
#> Call:
#> survival::coxph(formula = surv_obj ~ pdata[[variable]], data = pdata)
#> 
#>   n= 367, number of events= 144 
#> 
#>                        coef exp(coef)  se(coef)      z Pr(>|z|)   
#> pdata[[variable]] -0.018491  0.981679  0.007099 -2.605   0.0092 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                   exp(coef) exp(-coef) lower .95 upper .95
#> pdata[[variable]]    0.9817      1.019    0.9681    0.9954
#> 
#> Concordance= 0.569  (se = 0.026 )
#> Likelihood ratio test= 6.61  on 1 df,   p=0.01
#> Wald test            = 6.78  on 1 df,   p=0.009
#> Score (logrank) test = 6.81  on 1 df,   p=0.009
#> 
#> 
#> $summary_binary_variable
#>  Low High 
#>  257  110 
#> 
#> $cox_binary_object
#> Call:
#> survival::coxph(formula = surv_obj ~ pdata[[variable2]], data = pdata)
#> 
#>   n= 367, number of events= 144 
#> 
#>                           coef exp(coef) se(coef)      z Pr(>|z|)    
#> pdata[[variable2]]High -0.7056    0.4938   0.2058 -3.429 0.000605 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                        exp(coef) exp(-coef) lower .95 upper .95
#> pdata[[variable2]]High    0.4938      2.025    0.3299    0.7391
#> 
#> Concordance= 0.567  (se = 0.02 )
#> Likelihood ratio test= 13.32  on 1 df,   p=3e-04
#> Wald test            = 11.76  on 1 df,   p=6e-04
#> Score (logrank) test = 12.25  on 1 df,   p=5e-04
#> 
#> 
table(sig_stad2$TMEscore_CIR_binary)
#> 
#>  Low High 
#>  257  110 
# }
```
