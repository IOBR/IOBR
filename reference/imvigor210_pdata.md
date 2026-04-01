# IMvigor210 Bladder Cancer Immunotherapy Cohort Data

Clinical and biomarker data from the IMvigor210 clinical trial cohort.
Includes treatment response, survival outcomes, and immune biomarker
measurements for bladder cancer patients treated with atezolizumab.

## Usage

``` r
data(imvigor210_pdata)
```

## Format

A data frame with patients as rows and variables as columns:

- ID:

  Patient sample identifier

- BOR:

  Best overall response (CR, PR, SD, PD, NA)

- BOR_binary:

  Binary response classification (R=responder, NR=non-responder)

- OS_days:

  Overall survival time in days

- OS_status:

  Overall survival status (0=alive, 1=dead)

- Mutation_Load:

  Tumor mutation burden

- Neo_antigen_Load:

  Neoantigen load

- CD_8_T_effector:

  CD8+ T effector signature score

- Immune_Checkpoint:

  Immune checkpoint signature score

- Pan_F_TBRs:

  Pan-fibroblast TGF-β response signature

- Mismatch_Repair:

  Mismatch repair status or signature

- TumorPurity:

  Estimated tumor purity

## Source

IMvigor210 clinical trial (NCT02108652)

## References

Mariathasan S et al. TGFβ attenuates tumour response to PD-L1 blockade
by contributing to exclusion of T cells. Nature 554, 544-548 (2018).
doi:10.1038/nature25501
