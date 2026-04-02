# NULL Model Coefficients for MCPcounter

NULL Model Coefficients for MCPcounter

## Usage

``` r
data(null_models)
```

## Format

A \`data.frame\` with cell types in rows and coefficients in columns.

## Examples

``` r
data(null_models)
head(null_models)
#>                        Cell.population mu.133P2 sigma.133P2  mu.133A sigma.133A
#> T cells                        T cells 5.029214   0.4773043 6.598297  0.2664932
#> CD8 T cells                CD8 T cells 5.066625   0.4951542 6.715791  0.2184927
#> Cytotoxic lymphocytes        Cytotoxic 5.087922   0.7909333 6.267469  0.5181301
#> NK cells                      NK cells 4.882433   0.4123479 6.497035  0.3063453
#> B lineage                    B_derived 5.019908   0.5434946 6.452053  0.4779427
#> Monocytic lineage     Monocyte_derived 5.319001   0.4078621 7.022752  0.2821686
#>                         mu.HG1 sigma.HG1
#> T cells               5.614915 0.4028037
#> CD8 T cells           5.838888 0.4003190
#> Cytotoxic lymphocytes 5.013255 0.4994423
#> NK cells              5.473830 0.4332493
#> B lineage             5.418311 0.4228024
#> Monocytic lineage     5.236936 0.7015838
```
