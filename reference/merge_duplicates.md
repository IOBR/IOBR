# Merging the duplicates from the input matrix.

In case there are some duplicate rownames in the input matrix (*mat*),
this function will return a similar matrix with unique rownames and the
duplicate cases will be merged together based on the median values. When
*warn* is true a warning message will be written in case there are
duplicates, this warning message will include the *in_type* string to
indicate the current type of the input matrix.

## Usage

``` r
merge_duplicates(mat, warn = TRUE, in_type = NULL)
```

## Arguments

- mat:

  matrix

- warn:

  warning

- in_type:

  default is null
