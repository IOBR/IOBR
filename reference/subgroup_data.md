# Example Dataset for Subgroup Survival Analysis

An example dataset demonstrating the data structure required for
subgroup survival analysis using the
[`subgroup_survival`](https://iobr.github.io/IOBR/reference/subgroup_survival.md)
function. Contains simulated clinical and biomarker data with survival
outcomes.

## Usage

``` r
data(subgroup_data)
```

## Format

A data frame with clinical variables and biomarker scores:

- Patient_ID:

  Unique patient identifier

- ProjectID:

  Study or dataset identifier (e.g., "Dataset1")

- AJCC_stage:

  AJCC pathological stage (2, 3, 4)

- status:

  Event status (0=censored, 1=event)

- time:

  Follow-up time in months

- score:

  Continuous biomarker score (numeric values)

- score_binary:

  Binary biomarker classification ("High", "Low")

## Examples

``` r
data(subgroup_data)
head(subgroup_data)
#>   Patient_ID ProjectID AJCC_stage status   time       score score_binary
#> 1          1  Dataset1          2      0  88.73  0.60585688          Low
#> 2          2  Dataset1          2      0  88.23  0.73717229         High
#> 3          3  Dataset1          2      0  88.23 -0.35452887          Low
#> 4          4  Dataset1          2      0 105.70  0.79880007         High
#> 5          5  Dataset1          3      0 105.53 -0.09554256          Low
#> 6          6  Dataset1          2      1  25.50 -0.51527214          Low
```
