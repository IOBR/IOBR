# General RNA-seq Annotation

Generic gene annotation for RNA-seq data analysis. Provides basic gene
identifier mappings for various RNA-seq analysis workflows.

## Usage

``` r
data(anno_rnaseq)
```

## Format

A data frame with columns:

- probe_id:

  Gene identifier (platform specific).

- symbol:

  Official gene symbol.

## Source

Curated from multiple public RNA-seq resources

## Examples

``` r
data(anno_rnaseq)
head(anno_rnaseq)
#> # A tibble: 6 × 2
#>   probe_id        symbol      
#>   <fct>           <fct>       
#> 1 ENSG00000223972 DDX11L1     
#> 2 ENSG00000227232 WASH7P      
#> 3 ENSG00000278267 MIR6859-3   
#> 4 ENSG00000243485 RP11-34P13.3
#> 5 ENSG00000274890 MIR1302-9   
#> 6 ENSG00000237613 FAM138A     
```
