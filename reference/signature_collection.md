# Gene signature collection for pathway and immune analysis

A named list of gene signatures used in the IOBR package for immune
deconvolution, pathway scoring, functional annotation, and tumour
microenvironment (TME) characterization. Each element corresponds to a
predefined biological signature and contains a character vector of HGNC
gene symbols.

## Usage

``` r
data(signature_collection)
```

## Format

A named list of length 323. Each element is a character vector of gene
symbols. Representative entries include:

- CD_8_T_effector:

  Markers of CD8\\^{+}\\ effector T cells.

- DDR:

  DNA damage response and repair genes.

- Immune_Checkpoint:

  Immune checkpoint molecules.

- CellCycle_Reg:

  Core regulators of cell-cycle progression.

- Mismatch_Repair:

  Mismatch-repair pathway genes.

- TMEsocreA_CIR:

  TME-related signature used in TMEscore analysis.

- ...:

  Additional signatures are included in the list but are not
  individually listed here; all follow the same structure.
