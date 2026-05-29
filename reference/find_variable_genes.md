# Identify Variable Genes in Expression Data

Identifies variable genes from a gene expression dataset using specified
selection criteria. Supports multiple methods, including expression
thresholding and variability estimation via median absolute deviation
(MAD).

## Usage

``` r
find_variable_genes(
  eset,
  data_type = c("count", "normalized"),
  methods = c("low", "mad"),
  prop = 0.7,
  quantile = 0.75,
  min.mad = 0.1,
  feas = NULL
)
```

## Arguments

- eset:

  Numeric matrix. Gene expression data (genes as rows, samples as
  columns).

- data_type:

  Character. Type of data: \`"count"\` or \`"normalized"\`. Default is
  \`"count"\`.

- methods:

  Character vector. Methods for gene selection: \`"low"\`, \`"mad"\`.
  Default is \`c("low", "mad")\`.

- prop:

  Numeric. Proportion of samples in which a gene must be expressed.
  Default is 0.7.

- quantile:

  Numeric. Quantile threshold for minimum MAD (0.25, 0.5, 0.75). Default
  is 0.75.

- min.mad:

  Numeric. Minimum allowable MAD value. Default is 0.1.

- feas:

  Character vector or \`NULL\`. Additional features to include. Default
  is \`NULL\`.

## Value

Matrix subset of \`eset\` containing variable genes.

## Author

Dongqiang Zeng

## Examples

``` r
# Simulate data
set.seed(123)
sim_eset <- matrix(rnorm(100 * 20), 100, 20)
rownames(sim_eset) <- paste0("Gene", 1:100)
colnames(sim_eset) <- paste0("Sample", 1:20)

# Identify variable genes
eset_var <- find_variable_genes(
  eset = sim_eset,
  data_type = "normalized",
  methods = "mad",
  quantile = 0.25
)
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ min.mad = 0.1
#> ℹ Range of MAD: 0.41 to 1.72
#> ℹ 25% of variables will be filtered out...
head(eset_var)
#>           Sample1     Sample2    Sample3    Sample4     Sample5     Sample6
#> Gene3  1.55870831 -0.24669188 -0.2651451 -0.9385387 -0.63474826  1.02678506
#> Gene4  0.07050839 -0.34754260  0.5431941 -1.0525133 -0.02884155  0.75106130
#> Gene5  0.12928774 -0.95161857 -0.4143399 -0.4371595  0.67069597 -1.50916654
#> Gene6  1.71506499 -0.04502772 -0.4762469  0.3311792 -1.65054654 -0.09514745
#> Gene7  0.46091621 -0.78490447 -0.7886028 -2.0142105 -0.34975424 -0.89594782
#> Gene8 -1.26506123 -1.66794194 -0.5946173  0.2119804  0.75640644 -2.07075107
#>           Sample7    Sample8     Sample9    Sample10    Sample11    Sample12
#> Gene3 -0.03333034 -0.6930946  0.85520221  0.29959368 -0.01798024 -0.93656903
#> Gene4 -1.51606762  0.1188494  1.15293623  1.63905191 -0.13217513 -1.40078743
#> Gene5  0.79038534 -1.3647095  0.27627456  1.08461701 -2.54934277  0.16027754
#> Gene6 -0.21073418  0.5899827  0.14410466 -0.62456747  1.04057346 -0.27396237
#> Gene7 -0.65674293  0.2893440 -0.07562508  0.82592290  0.24972574 -0.98553911
#> Gene8 -1.41202579 -0.9042150  2.16141585 -0.04856835  2.41620737  0.08393068
#>         Sample13   Sample14   Sample15   Sample16   Sample17   Sample18
#> Gene3  0.8515247 -1.1477708  0.7148484 -0.9020980 -0.4539977 -3.0478609
#> Gene4 -0.7479300  0.3543522 -0.4316611  0.6270687 -0.5938646  1.8686555
#> Gene5  0.6302398  0.4247998  0.2276149  1.1203550 -1.7103797  1.7904242
#> Gene6  1.0966616  0.6483474  1.2949458  2.1272136 -0.2094484 -1.1010817
#> Gene7 -0.9884429 -1.2198100  0.5783349  0.3661144  2.4787458 -0.1681075
#> Gene8  1.1079950  0.1072350  1.3646728 -0.8747814  0.9897022  1.3752753
#>         Sample19    Sample20
#> Gene3 -0.0573241 -0.03265845
#> Gene4  1.2567478  1.63675735
#> Gene5  1.5874541 -0.32904197
#> Gene6  0.3194815 -2.60403817
#> Gene7  0.3815916  0.51398379
#> Gene8 -0.2436449 -0.88646801
```
