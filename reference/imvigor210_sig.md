# IMvigor210 Bladder Cancer Cohort Multi-omics Signatures

Comprehensive multi-omics dataset from the IMvigor210 phase II clinical
trial of metastatic urothelial cancer patients treated with atezolizumab
(anti-PD-L1). This dataset includes CIBERSORT immune cell deconvolution
scores, gene expression values, and various immune-related molecular
signatures.

## Usage

``` r
data(imvigor210_sig)
```

## Format

A data frame with 348 rows (patients) and 456 variables:

- ID – unique sample identifier ("SAM\*\*\*")

- 22 CIBERSORT immune-cell proportions (0–1 scale)

- immune-gene expression values

- immune signature scores (GEP, cytolytic activity, etc.)

- clinical / molecular features extracted from the original publication

## Source

IMvigor210 clinical trial (NCT02108652) Supplementary data from:
Mariathasan S et al. Nature 554, 544-548 (2018)
<https://www.nature.com/articles/nature25501>
<https://clinicaltrials.gov/ct2/show/NCT02108652>

## References

Mariathasan S et al. TGFβ attenuates tumour response to PD-L1 blockade
by contributing to exclusion of T cells. Nature 554, 544-548 (2018).
doi:10.1038/nature25501

Rosenberg JE et al. Atezolizumab in patients with locally advanced and
metastatic urothelial carcinoma who have progressed following treatment
with platinum-based chemotherapy: a single-arm, multicentre, phase 2
trial. Lancet 387, 1909-1920 (2016). doi:10.1016/S0140-6736(16)00561-4

Newman AM et al. Robust enumeration of cell subsets from tissue
expression profiles. Nature Methods 12, 453-457 (2015).
doi:10.1038/nmeth.3337

## Examples

``` r
data(imvigor210_sig)
head(imvigor210_sig)
#> # A tibble: 6 × 456
#>   ID        B_cells_naive_CIBERS…¹ B_cells_memory_CIBER…² Plasma_cells_CIBERSORT
#>   <chr>                      <dbl>                  <dbl>                  <dbl>
#> 1 SAMf2ce1…                 0.0290                0.0457                 0      
#> 2 SAM698d8…                 0.0812                0.00127                0      
#> 3 SAMc1b27…                 0.0124                0                      0.00127
#> 4 SAM85e41…                 0                     0                      0.00119
#> 5 SAMf275e…                 0                     0                      0.00949
#> 6 SAM7f0d9…                 0.0582                0.129                  0      
#> # ℹ abbreviated names: ¹​B_cells_naive_CIBERSORT, ²​B_cells_memory_CIBERSORT
#> # ℹ 452 more variables: T_cells_CD8_CIBERSORT <dbl>,
#> #   T_cells_CD4_naive_CIBERSORT <dbl>,
#> #   T_cells_CD4_memory_resting_CIBERSORT <dbl>,
#> #   T_cells_CD4_memory_activated_CIBERSORT <dbl>,
#> #   T_cells_follicular_helper_CIBERSORT <dbl>,
#> #   `T_cells_regulatory_(Tregs)_CIBERSORT` <dbl>, …
```
