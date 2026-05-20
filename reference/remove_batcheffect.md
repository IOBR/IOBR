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
  id_type,
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
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
eset_blca <- load_data("eset_blca")
#> ℹ Trying mirror 1/11: <https://github.com>
#> ✔ Download complete: "eset_blca"
# \donttest{
eset_corrected <- remove_batcheffect(
  eset_stad[1:1000, 1:5], eset_blca[1:1000, 1:5],
  id_type = "ensembl",
  data_type = "count"
)
#> ℹ The two expression matrices share 995 features
#> ℹ Processing method: sva::ComBat_seq
#> Found 2 batches
#> Using null model in ComBat-seq.
#> Adjusting for 0 covariate(s) or covariate level(s)
#> Estimating dispersions
#> Fitting the GLM model
#> Shrinkage off - using GLM estimates for parameters
#> Adjusting the data
#> ℹ Loading cached data: "anno_grch38"
#> ℹ Using local annotation (anno_grch38) for TPM conversion
#> ℹ No duplicate gene symbols found.
#> ✔ Applied log2 transformation
#> ℹ Count data processed with ComBat_seq
#> ✔ Applied log2 transformation
#> 
#> eset1 eset2 
#>     5     5 
#> >>== colors for group: 
#> >>== #374E55FF>>== #DF8F44FF
#> ✔ Applied log2 transformation
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
# }
```
