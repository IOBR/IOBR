# Deconvolve Immune Microenvironment Using xCell

Estimates immune cell fractions using the xCell algorithm. xCell
provides cell type enrichment scores for 64 immune and stromal cell
types from gene expression data.

## Usage

``` r
deconvo_xcell(eset, project = NULL, arrays = FALSE)
```

## Arguments

- eset:

  Gene expression matrix with HGNC gene symbols as row names and samples
  as columns.

- project:

  Optional project name to add as \`ProjectID\` column. Default is
  \`NULL\`.

- arrays:

  Logical indicating microarray data (\`TRUE\`) or RNA-seq (\`FALSE\`).
  Default is \`FALSE\`.

## Value

Data frame with xCell enrichment scores. Cell type columns are suffixed
with \`\_xCell\`.

## Author

Dongqiang Zeng

## Examples

``` r
xcell <- load_data("xCell.data")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "xCell.data"
if (!is.null(xcell)) {
  set.seed(123)
  sim_eset <- matrix(rnorm(length(xcell$genes) * 3), length(xcell$genes), 3)
  rownames(sim_eset) <- xcell$genes
  colnames(sim_eset) <- paste0("Sample", 1:3)
  
  # Run deconvolution
  result <- deconvo_xcell(eset = sim_eset, project = "TCGA-STAD")
  if (!is.null(result)) head(result)
}
#> ℹ Running xCell deconvolution
#> ℹ Loading cached data: "xCell.data"
#> ℹ Number of genes: 10808
#> ℹ GSVA version 2.6.2
#> ℹ Searching for rows with constant values
#> ℹ Calculating GSVA ranks
#> ℹ kcdf='auto' (default)
#> ℹ GSVA dense (classical) algorithm
#> ℹ Row-wise ECDF estimation with Gaussian kernels
#> ℹ Calculating row ECDFs
#> ℹ Calculating column ranks
#> ℹ GSVA dense (classical) algorithm
#> ℹ Calculating GSVA scores for 489 gene sets
#> ✔ Calculations finished
#>        ID ProjectID Adipocytes_xCell Astrocytes_xCell B-cells_xCell
#> 1 Sample1 TCGA-STAD     5.839975e-24     1.357894e-23  2.908985e-09
#> 2 Sample2 TCGA-STAD     8.783395e-23     0.000000e+00  0.000000e+00
#> 3 Sample3 TCGA-STAD     1.597735e-22     9.588740e-22  2.106454e-23
#>   Basophils_xCell CD4+_T-cells_xCell CD4+_Tcm_xCell CD4+_Tem_xCell
#> 1    9.765242e-07       0.000000e+00    6.85281e-24   7.074759e-08
#> 2    0.000000e+00       0.000000e+00    4.28731e-23   1.174048e-06
#> 3    1.886662e-06       1.708436e-22    0.00000e+00   0.000000e+00
#>   CD4+_memory_T-cells_xCell CD4+_naive_T-cells_xCell CD8+_T-cells_xCell
#> 1              0.000000e+00             1.323389e-22       0.000000e+00
#> 2              7.039736e-11             0.000000e+00       1.240735e-23
#> 3              2.904115e-23             1.363629e-22       0.000000e+00
#>   CD8+_Tcm_xCell CD8+_Tem_xCell CD8+_naive_T-cells_xCell    CLP_xCell
#> 1   0.000000e+00              0             1.691624e-23 4.020083e-23
#> 2   1.239081e-08              0             8.154746e-09 3.718837e-26
#> 3   4.224813e-23              0             8.740207e-23 7.997449e-23
#>      CMP_xCell Chondrocytes_xCell Class-switched_memory_B-cells_xCell
#> 1 0.000000e+00       5.717669e-24                        1.045045e-23
#> 2 4.917169e-10       3.528086e-22                        4.870414e-08
#> 3 0.000000e+00       0.000000e+00                        0.000000e+00
#>       DC_xCell Endothelial_cells_xCell Eosinophils_xCell Epithelial_cells_xCell
#> 1 0.000000e+00            1.760751e-23      1.935509e-23           0.000000e+00
#> 2 5.155503e-23            0.000000e+00      0.000000e+00           8.454114e-24
#> 3 1.051736e-23            0.000000e+00      0.000000e+00           3.173240e-22
#>   Erythrocytes_xCell Fibroblasts_xCell    GMP_xCell    HSC_xCell
#> 1       5.462807e-23      7.426554e-07 4.347742e-24 2.320835e-08
#> 2       0.000000e+00      2.609799e-21 0.000000e+00 0.000000e+00
#> 3       2.883233e-23      5.902750e-22 0.000000e+00 1.455165e-23
#>   Hepatocytes_xCell Keratinocytes_xCell    MEP_xCell    MPP_xCell    MSC_xCell
#> 1      3.375462e-08        1.723631e-24 1.022884e-06 0.000000e+00 5.958079e-23
#> 2      0.000000e+00        5.983594e-25 1.881100e-26 3.479724e-24 5.772338e-22
#> 3      0.000000e+00        1.258025e-23 1.336678e-07 0.000000e+00 1.331732e-22
#>   Macrophages_xCell Macrophages_M1_xCell Macrophages_M2_xCell Mast_cells_xCell
#> 1                 0         3.442342e-24         0.000000e+00     7.451198e-09
#> 2                 0         4.466303e-23         7.547952e-26     0.000000e+00
#> 3                 0         0.000000e+00         2.722865e-25     0.000000e+00
#>   Megakaryocytes_xCell Melanocytes_xCell Memory_B-cells_xCell
#> 1         0.000000e+00                 0          5.83088e-24
#> 2         0.000000e+00                 0          0.00000e+00
#> 3         3.501945e-24                 0          0.00000e+00
#>   Mesangial_cells_xCell Monocytes_xCell Myocytes_xCell NK_cells_xCell
#> 1          1.495642e-23    1.423937e-22   3.757921e-26   9.290121e-23
#> 2          3.803919e-22    0.000000e+00   1.841161e-23   0.000000e+00
#> 3          3.967947e-22    0.000000e+00   0.000000e+00   2.168879e-22
#>      NKT_xCell Neurons_xCell Neutrophils_xCell Osteoblast_xCell Pericytes_xCell
#> 1 3.203153e-23             0      0.000000e+00     1.215375e-22    9.117077e-07
#> 2 0.000000e+00             0      4.362301e-22     4.316634e-06    2.156544e-21
#> 3 1.080383e-22             0      0.000000e+00     8.312306e-07    1.613775e-21
#>   Plasma_cells_xCell Platelets_xCell Preadipocytes_xCell Sebocytes_xCell
#> 1       0.000000e+00    0.000000e+00        6.436607e-07    0.000000e+00
#> 2       5.369871e-08    1.706983e-22        0.000000e+00    0.000000e+00
#> 3       6.201039e-25    6.337865e-23        8.211132e-22    4.585329e-24
#>   Skeletal_muscle_xCell Smooth_muscle_xCell Tgd_cells_xCell Th1_cells_xCell
#> 1          3.006928e-23        4.818641e-23    0.000000e+00    0.000000e+00
#> 2          4.658381e-22        3.589720e-05    7.939448e-25    3.269045e-24
#> 3          0.000000e+00        1.985315e-05    0.000000e+00    4.155462e-21
#>   Th2_cells_xCell Tregs_xCell    aDC_xCell cDC_xCell   iDC_xCell
#> 1    3.066228e-24           0 1.701780e-10         0 1.04186e-08
#> 2    0.000000e+00           0 0.000000e+00         0 0.00000e+00
#> 3    0.000000e+00           0 3.699777e-22         0 0.00000e+00
#>   ly_Endothelial_cells_xCell mv_Endothelial_cells_xCell naive_B-cells_xCell
#> 1               0.000000e+00                          0        1.801548e-27
#> 2               0.000000e+00                          0        0.000000e+00
#> 3               3.806213e-23                          0        0.000000e+00
#>   pDC_xCell pro_B-cells_xCell ImmuneScore_xCell StromaScore_xCell
#> 1         0      1.517253e-26      6.906789e-09      3.713277e-07
#> 2         0      0.000000e+00      3.334617e-22      1.348816e-21
#> 3         0      4.354660e-27      2.795423e-22      3.750242e-22
#>   MicroenvironmentScore_xCell
#> 1                3.782345e-07
#> 2                1.682278e-21
#> 3                6.545665e-22
```
