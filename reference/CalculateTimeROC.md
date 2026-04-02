# Calculate Time-Dependent ROC Curve

Computes time-dependent ROC curve for survival models using the
\`timeROC\` package. Evaluates predictive accuracy at a specified time
quantile.

## Usage

``` r
CalculateTimeROC(model, newx, s, acture.y, modelname, time_prob = 0.9)
```

## Arguments

- model:

  A fitted survival model object.

- newx:

  A matrix or data frame of new data for prediction.

- s:

  Lambda value for prediction.

- acture.y:

  Data frame with \`time\` and \`status\` columns.

- modelname:

  Character string for model identification.

- time_prob:

  Numeric quantile for ROC calculation. Default is \`0.9\`.

## Value

An object of class \`timeROC\` containing ROC curve information.

## Author

Dongqiang Zeng

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE) &&
  requireNamespace("survival", quietly = TRUE) &&
  requireNamespace("timeROC", quietly = TRUE)) {
  library(survival)
  dat <- na.omit(lung[, c("time", "status", "age", "sex", "ph.ecog")])
  dat$status <- dat$status - 1
  x <- as.matrix(dat[, c("age", "sex", "ph.ecog")])
  y <- Surv(dat$time, dat$status)
  fit <- glmnet::glmnet(x, y, family = "cox")
  actual_outcome <- data.frame(time = dat$time, status = dat$status)
  roc_info <- CalculateTimeROC(
    model = fit, newx = x, s = 0.01, acture.y = actual_outcome,
    modelname = "glmnet Cox Model", time_prob = 0.5
  )
  print(roc_info$AUC)
}
#>       t=0     t=259 
#>        NA 0.6751424 
```
