# Reference profiles for B cell–related deconvolution (EPIC/IOBR)

`BRef` provides the B-cell–related reference gene expression profiles
used by the EPIC deconvolution model. Within the IOBR workflow, this
dataset is used to estimate immune cell fractions in bulk RNA-seq
samples. The dataset contains normalized gene expression profiles,
expression variability, and the signature genes required by EPIC.

## Usage

``` r
data(BRef)
```

## Format

A list of length 3:

- `refProfiles`: double matrix, approximately `49902 × 6`

- `refProfiles.var`: double matrix, approximately `49902 × 6`

- `sigGenes`: character vector, approximately `65` genes

## Details

This reference dataset includes the following components:

- `refProfiles`:

  A numeric matrix of size `nGenes × nCellTypes`. Rows correspond to
  genes, and columns correspond to B-cell–related subtypes. The values
  generally represent normalized expression (e.g., TPM). EPIC uses this
  matrix as the baseline expression profile for the reference cell types
  during deconvolution.

- `refProfiles.var`:

  A numeric matrix with identical dimensions as `refProfiles`. It
  represents gene-level expression variability for each cell type. These
  variability estimates are used as weights in EPIC’s weighted
  least-squares optimization. If variability is not used in a given
  workflow, EPIC assumes the same variability for all genes.

- `sigGenes`:

  A character vector listing the signature genes used by EPIC for B-cell
  deconvolution. Only these genes are included in the fitting procedure
  to improve the robustness and biological specificity of the model.
