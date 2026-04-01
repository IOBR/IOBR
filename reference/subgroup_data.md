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
