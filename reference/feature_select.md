# Feature Selection via Correlation or Differential Expression

Selects informative features using either correlation with a
quantitative response or differential expression (limma) for
binary/continuous responses.

## Usage

``` r
feature_select(
  x,
  y,
  method = c("cor", "dif"),
  family = c("spearman", "pearson"),
  cutoff = NULL,
  padjcut = NULL
)
```

## Arguments

- x:

  Numeric matrix. Features (rows) by samples (columns).

- y:

  Numeric or factor. Response vector (quantitative or binary).

- method:

  Character. "cor" (correlation) or "dif" (differential expression).
  Default c("cor","dif").

- family:

  Character. Correlation method if method = "cor": "spearman" or
  "pearson".

- cutoff:

  Numeric. Absolute correlation (for cor) or \|log2FC\| (for dif)
  threshold.

- padjcut:

  Numeric. Adjusted p-value cutoff.

## Value

Character vector of selected feature names.

## Examples

``` r
# Simulate data
set.seed(123)
sim_eset <- matrix(rnorm(100 * 50), 100, 50)
rownames(sim_eset) <- c("PDCD1", paste0("Gene", 2:100))
colnames(sim_eset) <- paste0("Sample", 1:50)
pd1 <- as.numeric(sim_eset["PDCD1", ])
group <- ifelse(pd1 > mean(pd1), "high", "low")

# Correlation method
pd1_cor <- feature_select(
  x = sim_eset, y = pd1, method = "cor",
  family = "pearson", padjcut = 0.05, cutoff = 0.5
)
#> Deteching more than two levels in y, we will treat y as a quantitative varibale

# Differential expression method
pd1_dif <- feature_select(
  x = sim_eset, y = pd1, method = "dif",
  padjcut = 0.05, cutoff = 2
)
#> Deteching more than two levels in y, we will treat y as a quantitative varibale
#> For quantitative varibale, upper 25% and bottom 25% samples
#>               were treated as upregulated group and downregulated group.
pd1_dif_2 <- feature_select(
  x = sim_eset, y = group,
  method = "dif", padjcut = 0.05, cutoff = 2
)
#> Deteching two levels in y, we will treat y as a binary varibale
```
