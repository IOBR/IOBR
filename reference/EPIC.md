# Estimate the proportion of immune and cancer cells.

`EPIC` takes as input bulk gene expression data (RNA-seq) and returns
the proportion of mRNA and cells composing the various samples.

## Usage

``` r
EPIC(
  bulk,
  reference = NULL,
  mRNA_cell = NULL,
  mRNA_cell_sub = NULL,
  sigGenes = NULL,
  scaleExprs = TRUE,
  withOtherCells = TRUE,
  constrainedSum = TRUE,
  rangeBasedOptim = FALSE
)
```

## Arguments

- bulk:

  A matrix (`nGenes` x `nSamples`) of the genes expression from each
  bulk sample (the counts should be given in TPM, RPKM or FPKM when
  using the prebuilt reference profiles). This matrix needs to have
  rownames telling the gene names (corresponds to the gene symbol in the
  prebuilt reference profiles (e.g. CD8A, MS4A1) - no conversion of IDs
  is performed at the moment). It is advised to keep all genes in the
  bulk instead of a subset of signature genes (except if
  `scaleExprs = FALSE` in which case it doesn't make any difference).

- reference:

  (optional): A string or a list defining the reference cells. It can
  take multiple formats: - \`NULL\`: to use the default reference
  profiles and genes signature `TRef` - A character: either `"BRef"` or
  `"TRef"` to use the reference cells and genes signature of the
  corresponding datasets (see `BRef` and `TRef`) - A list containing: -
  \`\$refProfiles\`: a matrix (`nGenes` x `nCellTypes`) of the reference
  cells genes expression (don't include a column of the 'other cells'
  (representing usually the cancer cells for which such a profile is
  usually not conserved between samples); the rownames needs to be
  defined as well as the colnames giving the names of each gene and
  reference cell types respectively. It is advised to keep all genes in
  this `refProfiles` matrix instead of a subset of signature genes -
  \`\$sigGenes\`: a character vector of the gene names to use as
  signature - sigGenes can also be given as a direct input to EPIC
  function - \`\$refProfiles.var\` (optional): a matrix (`nGenes` x
  `nCellTypes`) of the variability of each gene expression for each cell
  type, which is used to define weights on each gene for the
  optimization (if this is absent, we assume an identical variability
  for all genes in all cells) - it needs to have the same dimnames than
  refProfiles

- mRNA_cell:

  (optional): A named numeric vector: tells (in arbitrary units) the
  amount of mRNA for each of the reference cells and of the other
  uncharacterized (cancer) cell. Two names are of special meaning:
  *"otherCells"* - used for the mRNA/cell value of the "other cells"
  from the sample (i.e. the cell type that don't have any reference gene
  expression profile) ; and *default* - used for the mRNA/cell of the
  cells from the reference profiles for which no specific value is given
  in mRNA_cell (i.e. if mRNA_cell=c(Bcells=2, NKcells=2.1,
  otherCells=3.5, default=1), then if the refProfiles described Bcells,
  NKcells and Tcells, we would use a value of 3.5 for the "otherCells"
  that didn't have any reference profile and a default value of 1 for
  the Tcells when computing the cell fractions). To note: if data is in
  tpm, this mRNA per cell would ideally correspond to some number of
  transcripts per cell.

- mRNA_cell_sub:

  (optional): This can be given instead of `mRNA_cell` (or in addition
  to it). It is also a named numeric vector, used to replace only the
  mRNA/cell values from some cell types (or to add values for new cell
  types). The values given in mRNA_cell_sub will overwrite the default
  values as well as those that might have been given by mRNA_cell.

- sigGenes:

  (optional): a character vector of the gene names to use as signature
  for the deconvolution. In principle this is given with the reference
  as the "reference\$sigGenes" but if we give a value for this input
  variable, it is these signature genes that will be used instead of the
  ones given with the reference profile.

- scaleExprs:

  (optional, default is TRUE): boolean telling if the bulk samples and
  reference gene expression profiles should be rescaled based on the
  list of genes in common between the them (such a rescaling is
  recommanded).

- withOtherCells:

  (optional, default is TRUE): if EPIC should allow for an additional
  cell type for which no gene expression reference profile is available
  or if the bulk is assumed to be composed only of the cells with
  reference profiles.

- constrainedSum:

  (optional, default is TRUE): tells if the sum of all cell types should
  be constrained to be \< 1. When `withOtherCells=FALSE`, there is
  additionally a constrain the the sum of all cell types with reference
  profiles must be \> 0.99.

- rangeBasedOptim:

  (optional): when this is FALSE (the default), the least square
  optimization is performed as described in [Racle et al., 2017,
  eLife](https://elifesciences.org/articles/26476), which is
  recommanded. When this variable is TRUE, EPIC uses the variability of
  each gene from the reference profiles in another way: instead of
  defining weights (based on the variability) for the fit of each gene,
  we define a range of values accessible for each gene (based on the
  gene expression value in the reference profile +/- the variability
  values). The error that the optimization tries to minimize is by how
  much the predicted gene expression is outside of this allowed range of
  values.

## Value

A list of 3 matrices: - \`mRNAProportions\`: (`nSamples` x
(`nCellTypes+1`)) the proportion of mRNA coming from all cell types with
a ref profile + the uncharacterized other cell - \`cellFractions\`:
(`nSamples` x (`nCellTypes+1`)) this gives the proportion of cells from
each cell type after accounting for the mRNA / cell value - \`fit.gof\`:
(`nSamples` x 12) a matrix telling the quality for the fit of the
signature genes in each sample. It tells if the minimization converged,
and other info about this fit comparing the measured gene expression in
the sigGenes vs predicted gene expression in the sigGenes

## Details

This function uses a constrained least square minimization to estimate
the proportion of each cell type with a reference profile and another
uncharacterized cell type in bulk gene expression samples.

The names of the genes in the bulk samples, the reference samples and in
the gene signature list need to be the same format (gene symbols are
used in the predefined reference profiles). The full list of gene names
don't need to be exactly the same between the reference and bulk
samples: *EPIC* will use the intersection of the genes. In case of
duplicate gene names, *EPIC* will use the median value per duplicate -
if you want to consider these cases differently, you can remove the
duplicates before calling *EPIC*.

## Examples

``` r
# \donttest{
melanoma_data <- load_data("melanoma_data")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "melanoma_data"
TRef <- load_data("TRef")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "TRef"
res1 <- EPIC(melanoma_data$counts)
#> ℹ Loading cached data: "TRef"
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "mRNA_cell_default"
#> Warning: The optimization didn't fully converge for some samples:
#> LAU1255; LAU355
#>  - check fit.gof for the convergeCode and convergeMessage
#> Warning: mRNA_cell value unknown for some cell types: CAFs, Endothelial - using the default value of 0.4 for these but this might bias the true cell proportions from all cell types.
res1$cellFractions
#>             Bcells         CAFs CD4_Tcells  CD8_Tcells  Endothelial Macrophages
#> LAU125  0.01016774 2.359222e-04 0.03041567 0.009302975 2.554933e-02  0.01193535
#> LAU1255 0.04136029 5.590262e-04 0.05849291 0.129023216 4.649865e-07  0.01961274
#> LAU1314 0.67562140 1.163239e-08 0.08198931 0.015668658 2.799781e-07  0.00144306
#> LAU355  0.45493081 7.457349e-09 0.26813661 0.016429308 4.095153e-09  0.00861752
#>              NKcells otherCells
#> LAU125  9.058697e-09  0.9123930
#> LAU1255 4.706589e-07  0.7509509
#> LAU1314 4.095899e-03  0.2211814
#> LAU355  4.601164e-08  0.2518857
res2 <- EPIC(melanoma_data$counts, TRef)
#> ℹ Loading cached data: "mRNA_cell_default"
#> Warning: The optimization didn't fully converge for some samples:
#> LAU1255; LAU355
#>  - check fit.gof for the convergeCode and convergeMessage
#> Warning: mRNA_cell value unknown for some cell types: CAFs, Endothelial - using the default value of 0.4 for these but this might bias the true cell proportions from all cell types.
res3 <- EPIC(bulk = melanoma_data$counts, reference = TRef)
#> ℹ Loading cached data: "mRNA_cell_default"
#> Warning: The optimization didn't fully converge for some samples:
#> LAU1255; LAU355
#>  - check fit.gof for the convergeCode and convergeMessage
#> Warning: mRNA_cell value unknown for some cell types: CAFs, Endothelial - using the default value of 0.4 for these but this might bias the true cell proportions from all cell types.
res4 <- EPIC(melanoma_data$counts, reference = "TRef")
#> ℹ Loading cached data: "TRef"
#> ℹ Loading cached data: "mRNA_cell_default"
#> Warning: The optimization didn't fully converge for some samples:
#> LAU1255; LAU355
#>  - check fit.gof for the convergeCode and convergeMessage
#> Warning: mRNA_cell value unknown for some cell types: CAFs, Endothelial - using the default value of 0.4 for these but this might bias the true cell proportions from all cell types.
res5 <- EPIC(melanoma_data$counts, mRNA_cell_sub = c(
  Bcells = 1,
  otherCells = 5
))
#> ℹ Loading cached data: "TRef"
#> ℹ Loading cached data: "mRNA_cell_default"
#> Warning: The optimization didn't fully converge for some samples:
#> LAU1255; LAU355
#>  - check fit.gof for the convergeCode and convergeMessage
#> Warning: mRNA_cell value unknown for some cell types: CAFs, Endothelial - using the default value of 0.4 for these but this might bias the true cell proportions from all cell types.
# Various possible ways of calling EPIC function. res 1 to 4 should
# give exactly the same outputs, and the elements res1$cellFractions
# should be equal to the example predictions found in
# melanoma_data$cellFractions.pred for these first 4 results.
# The values of cellFraction for res5 will be different due to the use of
# other mRNA per cell values for the B and other cells.
# }
```
