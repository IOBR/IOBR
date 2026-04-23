# Tumor Microenvironment (TME) Deconvolution Pipeline

Executes an integrated TME analysis on a gene expression matrix:
performs immune/stromal cell deconvolution using multiple algorithms,
computes signature scores, and aggregates results. Designed for
exploratory immunogenomic profiling.

## Usage

``` r
iobr_deconvo_pipeline(
  eset,
  project,
  array,
  tumor_type,
  path = NULL,
  permutation = 1000
)
```

## Arguments

- eset:

  Numeric matrix. Gene expression (TPM/log scale) with genes in rows.

- project:

  Character. Project name (used in output naming).

- array:

  Logical. Whether data originated from an array platform. Affects
  deconvolution choices.

- tumor_type:

  Character. Tumor type code (e.g., "stad") used by certain methods.

- path:

  Character. Output directory. Default is NULL (uses tempdir()).

- permutation:

  Integer. Number of permutations for CIBERSORT (and similar). Default
  is 1000.

## Value

Data frame integrating cell fractions and signature scores (also writes
intermediate outputs to disk).

## Author

Dongqiang Zeng

## Examples

``` r
eset_stad <- load_data("eset_stad")
#> ℹ Loading cached data: "eset_stad"
anno_grch38 <- load_data("anno_grch38")
#> ℹ Loading cached data: "anno_grch38"
eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#> ℹ Row number of original eset: 60483
#> ✔ 100% of probes in expression set were annotated
#> ℹ Found 2293 duplicate symbols, using "mean" method
#> ℹ Row number after filtering duplicated gene symbol: 50181
# \donttest{
res <- iobr_deconvo_pipeline(
  eset = eset, project = "STAD",
  array = FALSE, tumor_type = "stad",
  path = tempdir(), permutation = 10
)
#> ✔ Applied log2 transformation
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running CIBERSORT
#> ℹ Loading cached data: "lm22"
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running EPIC deconvolution
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "TRef"
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "mRNA_cell_default"
#> Warning: The optimization didn't fully converge for some samples:
#> TCGA-BR-6455; TCGA-FP-7916
#>  - check fit.gof for the convergeCode and convergeMessage
#> Warning: mRNA_cell value unknown for some cell types: CAFs, Endothelial - using the default value of 0.4 for these but this might bias the true cell proportions from all cell types.
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running MCP-counter deconvolution
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running xCell deconvolution
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "xCell.data"
#> ℹ Number of genes: 10783
#> ℹ GSVA version 2.4.9
#> ℹ Searching for rows with constant values
#> ℹ Calculating GSVA ranks
#> ℹ kcdf='auto' (default)
#> ℹ GSVA dense (classical) algorithm
#> ℹ Row-wise ECDF estimation with Gaussian kernels
#> ℹ Calculating row ECDFs
#> ℹ Calculating column ranks
#> ℹ GSVA dense (classical) algorithm
#> ℹ Calculating GSVA scores
#> ✔ Calculations finished
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running ESTIMATE
#> ℹ Loading cached data: "common_genes"
#> Merged dataset includes 10148 genes (264 mismatched).
#> ℹ Loading cached data: "SI_geneset"
#> 1 gene set: StromalSignature overlap=138
#> 2 gene set: ImmuneSignature overlap=140
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running TIMER deconvolution
#> ℹ Enter batch mode
#> ℹ Loading immune gene expression
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "immuneCuratedData"
#> ℹ Outlier genes: ACTB ACTG1 CD74 COL1A1 EEF1A1 ERBB2 FLNA IGHG1 IGKC MT-CO1 MT-CO2 MT-ND4 MT-RNR2 MYH11
#> ℹ Removing batch effects for stad
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "cancer_type_genes"
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running quanTIseq deconvolution
#> ℹ Running quanTIseq deconvolution module
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "quantiseq_data"
#> ℹ Gene expression normalization and re-annotation (arrays: FALSE)
#> ℹ Loading cached data: "quantiseq_data"
#> ℹ Removing 17 noisy genes
#> ℹ Removing 15 genes with high expression in tumors
#> ℹ Signature genes found in data set: 138/138 (100%)
#> ℹ Mixture deconvolution (method: lsei)
#> ✔ Deconvolution successful!
#> Warning: Data values appear small (< 50).
#> ℹ Input should be in TPM/FPKM scale, not log-transformed
#> ℹ Running IPS calculation
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "ips_gene_set"
#> [1] ">>>>> TME cell deconvolution was completed: STAD"
#> ℹ Calculating signature scores using PCA, z-score, and ssGSEA methods
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ Step 1/3: PCA method
#> ℹ Step 2/3: z-score method
#> ℹ Step 3/3: ssGSEA method
#> ℹ GSVA version 2.4.9
#> ℹ Searching for rows with constant values
#> ℹ Calculating ssGSEA scores for 280 gene sets
#> ℹ Calculating ranks
#> ℹ Calculating rank weights
#> ℹ Normalizing ssGSEA scores
#> ✔ Calculations finished
#> [1] ">>>>> Signature esitmation was completed: STAD"
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "hallmark"
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "go_bp"
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "go_cc"
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "go_mf"
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "kegg"
#> ℹ Trying mirror 1/4: <https://github.com>
#> ✔ Download complete: "reactome"
#> ℹ Calculating signature scores using ssGSEA method
#> ℹ Log2 transformation not necessary (data appears to already be log-scaled)
#> ℹ Calculating scores for 12059 signature(s)
#> ℹ GSVA version 2.4.9
#> ℹ Searching for rows with constant values
#> ℹ Calculating ssGSEA scores for 12059 gene sets
#> ℹ Calculating ranks
#> ℹ Calculating rank weights
#> Calculating ssGSEA scores ■■■■■■■■■■■■■■■■                  50% | ETA:  2s
#> Calculating ssGSEA scores ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s
#> ℹ Normalizing ssGSEA scores
#> ✔ Calculations finished
#> [1] ">>>>> HALLMARK GO KEGG REACTOME esitmation was completed: STAD"
# }
```
