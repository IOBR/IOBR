# Deconvolute Tumor Microenvironment Using TIMER

Performs deconvolution of the tumor microenvironment using the TIMER
algorithm. Processes multiple cancer datasets, removes batch effects,
and estimates immune cell type abundances.

## Usage

``` r
deconvolute_timer.default(args)
```

## Arguments

- args:

  List or environment containing parameters:

  outdir

  :   Character. Output directory path.

  batch

  :   Character. File containing paths to expression data and cancer
      types.

## Value

Matrix of abundance scores for different immune cell types across
multiple cancer samples.

## Examples

``` r
# \donttest{
# file
tf <- tempfile(fileext = ".csv")
write.table(data.frame("exp1", "luad", "exp2", "brca"),
  file = tf, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
)
outdir <- tempdir()
args <- list(outdir = outdir, batch = tf)
results <- deconvolute_timer.default(args)
#> ℹ Enter batch mode
#> ℹ Loading immune gene expression
#> ℹ Loading cached data: "immuneCuratedData"
#> Warning: cannot open file 'exp1': No such file or directory
#> Error in file(file, "rt"): cannot open the connection
# }
```
