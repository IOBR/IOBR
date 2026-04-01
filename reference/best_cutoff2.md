# Extract Best Cutoff and Add Binary Variable to Data Frame

Finds the best cutoff point for a continuous variable in survival
analysis. Takes input data containing a continuous variable and survival
information (time and status). Returns the modified input data with a
new binary variable created based on the best cutoff point.

## Usage

``` r
best_cutoff2(
  pdata,
  variable,
  time = "time",
  status = "status",
  print_result = TRUE
)
```

## Arguments

- pdata:

  Data frame containing survival information and features.

- variable:

  Character string specifying the continuous variable name.

- time:

  Character string specifying the time-to-event column name. Default is
  \`"time"\`.

- status:

  Character string specifying the event status column name. Default is
  \`"status"\`.

- print_result:

  Logical indicating whether to print results. Default is \`TRUE\`.

## Value

List containing:

- pdata:

  Data frame with binary variable added

- best_cutoff:

  Numeric cutoff value

- cox_continuous_object:

  Cox model summary for continuous variable

- summary_binary_variable:

  Summary of binary variable

- cox_binary_object:

  Cox model summary for binary variable

## Author

Dongqiang Zeng

## Examples

``` r
# \donttest{
sig_stad <- load_data("sig_stad")
result <- best_cutoff2(
  pdata = sig_stad,
  variable = "TMEscore_CIR",
  time = "OS_time",
  status = "OS_status",
  print_result = TRUE
)
#> ✔ Best cutoff for "TMEscore_CIR": 5.697
#> $best_cutoff
#> [1] "best cutoff = 5.69725525158544"
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
#>  256  111 
#> 
#> $cox_binary_object
#> Call:
#> survival::coxph(formula = surv_obj ~ pdata[[variable2]], data = pdata)
#> 
#>   n= 367, number of events= 144 
#> 
#>                           coef exp(coef) se(coef)     z Pr(>|z|)    
#> pdata[[variable2]]High -0.7057    0.4938   0.2058 -3.43 0.000604 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                        exp(coef) exp(-coef) lower .95 upper .95
#> pdata[[variable2]]High    0.4938      2.025    0.3299     0.739
#> 
#> Concordance= 0.567  (se = 0.02 )
#> Likelihood ratio test= 13.32  on 1 df,   p=3e-04
#> Wald test            = 11.76  on 1 df,   p=6e-04
#> Score (logrank) test = 12.25  on 1 df,   p=5e-04
#> 
#> 
result$best_cutoff
#> [1] 5.697255
# }
```
