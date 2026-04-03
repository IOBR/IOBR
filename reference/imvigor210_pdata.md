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

  Pan-fibroblast TGF-\\\beta\\ response signature

- Mismatch_Repair:

  Mismatch repair status or signature

- TumorPurity:

  Estimated tumor purity

## Source

IMvigor210 clinical trial (NCT02108652)

## References

Mariathasan S et al. TGF\\\beta\\ attenuates tumour response to PD-L1
blockade by contributing to exclusion of T cells. Nature 554, 544-548
(2018). doi:10.1038/nature25501

## Examples

``` r
data(imvigor210_pdata)
head(imvigor210_pdata)
#> # A tibble: 6 × 12
#>   ID           BOR   BOR_binary OS_days OS_status Mutation_Load Neo_antigen_Load
#>   <chr>        <chr> <chr>      <chr>   <chr>     <chr>         <chr>           
#> 1 SAM00b9e5c5… NA    NA         57.166… 1         NA            NA              
#> 2 SAM0257bbbb… SD    NR         469.15… 1         18            4.6862745099999…
#> 3 SAM025b45c2… PD    NR         263.16… 1         1             0.31372549      
#> 4 SAM032c6423… PD    NR         74.907… 1         44            6.1960784310000…
#> 5 SAM04c589eb… NA    NA         20.698… 0         50            NA              
#> 6 SAM0571f17f… SD    NR         136.01… 1         2             1.4705882349999…
#> # ℹ 5 more variables: CD_8_T_effector <dbl>, Immune_Checkpoint <dbl>,
#> #   Pan_F_TBRs <chr>, Mismatch_Repair <chr>, TumorPurity <dbl>
```
