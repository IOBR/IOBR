# Removing Batch Effect from Expression Sets

Removes batch effects from expression datasets using sva::ComBat (for
microarray/TPM data) or sva::ComBat_seq (for RNA-seq count data).
Generates PCA plots to visualize data before and after correction.

## Usage

``` r
remove_batcheffect(
  eset1,
  eset2,
  eset3 = NULL,
  id_type = "ensembl",
  data_type = c("array", "count", "tpm"),
  cols = "normal",
  palette = "jama",
  log2 = TRUE,
  check_eset = TRUE,
  adjust_eset = TRUE,
  repel = FALSE,
  path = NULL
)
```

## Arguments

- eset1:

  First expression set (matrix or data frame with genes as rows).

- eset2:

  Second expression set.

- eset3:

  Optional third expression set. Use \`NULL\` if not available.

- id_type:

  Type of gene ID in expression sets (e.g., \`"ensembl"\`,
  \`"symbol"\`). Required for count data normalization.

- data_type:

  Type of data: \`"array"\`, \`"count"\`, or \`"tpm"\`. Default is
  \`"array"\`.

- cols:

  Color scale for PCA plot. Default is \`"normal"\`.

- palette:

  Color palette for PCA plot. Default is \`"jama"\`.

- log2:

  Whether to perform log2 transformation. Default is \`TRUE\`. Ignored
  for count data.

- check_eset:

  Whether to check expression sets for errors. Default is \`TRUE\`.

- adjust_eset:

  Whether to adjust expression sets by removing problematic features.
  Default is \`TRUE\`.

- repel:

  Whether to add repelling labels to PCA plot. Default is \`FALSE\`.

- path:

  Directory where results should be saved. Default is \`NULL\` (display
  only).

## Value

Expression matrix after batch correction.

## References

Zhang Y, et al. ComBat-seq: batch effect adjustment for RNA-seq count
data. NAR Genomics and Bioinformatics. 2020;2(3):lqaa078.
doi:10.1093/nargab/lqaa078

Leek JT, et al. The sva package for removing batch effects and other
unwanted variation in high-throughput experiments. Bioinformatics.
2012;28(6):882-883.

## Author

Dongqiang Zeng

## Examples

``` r
# Simulate data
set.seed(123)
sim_eset1 <- matrix(rnorm(100 * 5, mean = 10, sd = 2), 100, 5)
sim_eset2 <- matrix(rnorm(100 * 5, mean = 12, sd = 2), 100, 5)
rownames(sim_eset1) <- rownames(sim_eset2) <- paste0("Gene", 1:100)
colnames(sim_eset1) <- paste0("S1_", 1:5)
colnames(sim_eset2) <- paste0("S2_", 1:5)

# Run batch correction
if (requireNamespace("sva", quietly = TRUE) && requireNamespace("BiocParallel", quietly = TRUE)) {
  eset_corrected <- remove_batcheffect(sim_eset1, sim_eset2, data_type = "tpm")
  if (!is.null(eset_corrected)) head(eset_corrected)
}
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ The two expression matrices share 100 features
#> ℹ Processing method: sva::ComBat
#> Found2batches
#> Adjusting for0covariate(s) or covariate level(s)
#> Standardizing Data across genes
#> Fitting L/S model and finding priors
#> Finding parametric adjustments
#> Adjusting the Data
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> 
#> eset1 eset2 
#>     5     5 
#> >>== colors for group: 
#> >>== #374E55FF>>== #DF8F44FF
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> 
#> eset1 eset2 
#>     5     5 
#> >>== colors for group: 
#> >>== #374E55FF>>== #DF8F44FF
#>            S1_1      S1_2     S1_3      S1_4      S1_5      S2_1      S2_2
#> Gene1  9.792976  9.944992 14.91214  9.895600 10.695508 10.109164 13.361392
#> Gene2  9.976222 11.449648 12.48123  9.522742  8.426746  9.566679 11.449648
#> Gene3 14.048781 11.021545 10.59517  9.683681  9.629662 12.978768 10.960286
#> Gene4 11.418995 10.960286 12.22496  9.456422 10.960286 12.224961  8.749894
#> Gene5 11.275189  9.350557 10.10916 10.238317 12.145525  8.260313 12.317743
#> Gene6 13.906184 11.275189 10.00584 11.545079  8.123103 10.730401 10.831304
#>            S2_3     S2_4      S2_5
#> Gene1  9.566679 11.54508  9.236308
#> Gene2  8.577297 10.10916  9.976222
#> Gene3  9.456422 12.14552 11.578457
#> Gene4 10.653537 12.59930 13.558965
#> Gene5  8.849318 11.36695 12.820878
#> Gene6 12.145525 11.27519  9.944992
```
