# Fix Expression Mixture Matrix

Processes expression matrix by mapping genes, converting log values,
quantile normalization (if arrays), and TPM normalization.

## Usage

``` r
fixMixture(mix.mat, arrays = FALSE)
```

## Arguments

- mix.mat:

  Expression matrix with genes as rows.

- arrays:

  Logical indicating if data is from arrays. Default is FALSE.

## Value

Processed expression matrix.
