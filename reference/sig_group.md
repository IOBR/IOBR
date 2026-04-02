# Grouped gene signatures for IOBR analysis

A named list that organizes gene signatures into functional or
biological categories. Each element of the list is a character vector
containing the names of gene signatures defined in
`signature_collection`. A total of 43 signature groups are included,
covering tumour intrinsic pathways, immune-related processes, stromal
activity, TME characteristics and immuno-oncology biomarkers. These
groups are used in IOBR to conveniently select sets of signatures for
scoring and visualization.

## Usage

``` r
data(sig_group)
```

## Format

A named list of length 43. Each element is a character vector of
signature names. Representative groups include:

- tumor_signature:

  Signatures related to intrinsic tumour biology such as cell cycle, DNA
  damage repair and histone regulation.

- EMT:

  Epithelial–mesenchymal transition (EMT)–associated signatures.

- io_biomarkers:

  Immuno-oncology biomarker–related signatures.

- immu_microenvironment:

  Immune microenvironment–related signatures.

- immu_suppression:

  Immune suppression–related signatures.

- immu_exclusion:

  Signatures associated with immune exclusion and stromal barriers.

- TCR_BCR:

  T-cell and B-cell receptor pathway signatures.

- tme_signatures1:

  Tumour microenvironment signature panel (set 1).

- tme_signatures2:

  Tumour microenvironment signature panel (set 2).

- ...:

  Additional groups are included (43 total), but not listed individually
  here; all groups follow the same structure.

## Examples

``` r
data(sig_group)
head(sig_group)
#> $tumor_signature
#>  [1] "CellCycle_Reg"                            
#>  [2] "Cell_cycle"                               
#>  [3] "DDR"                                      
#>  [4] "Mismatch_Repair"                          
#>  [5] "Histones"                                 
#>  [6] "Homologous_recombination"                 
#>  [7] "Nature_metabolism_Hypoxia"                
#>  [8] "Molecular_Cancer_m6A"                     
#>  [9] "MT_exosome"                               
#> [10] "Positive_regulation_of_exosomal_secretion"
#> [11] "Ferroptosis"                              
#> [12] "EV_Cell_2020"                             
#> 
#> $EMT
#> [1] "Pan_F_TBRs" "EMT1"       "EMT2"       "EMT3"       "WNT_target"
#> 
#> $io_biomarkers
#>  [1] "TMEscore_CIR"                    "TMEscoreA_CIR"                  
#>  [3] "TMEscoreB_CIR"                   "T_cell_inflamed_GEP_Ayers_et_al"
#>  [5] "CD_8_T_effector"                 "IPS_IPS"                        
#>  [7] "Immune_Checkpoint"               "Exhausted_CD8_Danaher_et_al"    
#>  [9] "Pan_F_TBRs"                      "Mismatch_Repair"                
#> [11] "APM"                            
#> 
#> $immu_microenvironment
#>  [1] "GPAGs"                       "PPAGs"                      
#>  [3] "Pan_F_TBRs"                  "MDSC_Peng_et_al"            
#>  [5] "TAM_Peng_et_al"              "MicroenvironmentScore_xCell"
#>  [7] "TMEscore"                    "ImmuneScore_xCell"          
#>  [9] "ImmuneScore_estimate"        "StromalScore_estimate"      
#> 
#> $immu_suppression
#> [1] "Pan_F_TBRs"                          
#> [2] "Fibroblasts_MCPcounter"              
#> [3] "Immune_Checkpoint"                   
#> [4] "Exhausted_CD8_Danaher_et_al"         
#> [5] "MDSC_Wang_et_al"                     
#> [6] "Macrophages_M2_cibersort"            
#> [7] "Tregs_quantiseq"                     
#> [8] "T_cells_regulatory_(Tregs)_CIBERSORT"
#> 
#> $immu_exclusion
#>  [1] "Pan_F_TBRs"                          
#>  [2] "CAF_Peng_et_al"                      
#>  [3] "Fibroblasts_MCPcounter"              
#>  [4] "CAFs_EPIC"                           
#>  [5] "EMT1"                                
#>  [6] "EMT2"                                
#>  [7] "WNT_target"                          
#>  [8] "TGFb_Family_Member_Receptor_Li_et_al"
#>  [9] "TGFb_Family_Member_Li_et_al"         
#> [10] "Macrophages_M2_CIBERSORT"            
#> [11] "Macrophages_M2_quantiseq"            
#> [12] "TAM_Peng_et_al"                      
#> 
```
