# Row Bind Multiple Data Sets

Combines two or three data frames or matrices vertically using
\`rbind\`. Ensures compatibility of input data before binding by
aligning columns.

## Usage

``` r
rbind_iobr(data1, data2, data3 = NULL)
```

## Arguments

- data1:

  Data frame or matrix. First dataset to combine.

- data2:

  Data frame or matrix. Second dataset to combine.

- data3:

  Data frame or matrix or \`NULL\`. Optional third dataset. Default is
  \`NULL\`.

## Value

Combined data frame resulting from row binding the input datasets.

## Author

Dongqiang Zeng

## Examples

``` r
data1 <- data.frame(A = 1:5, B = letters[1:5])
data2 <- data.frame(A = 6:10, B = letters[6:10])
combined_data <- rbind_iobr(data1, data2)

# With three datasets
data3 <- data.frame(A = 11:15, B = letters[11:15])
combined_data <- rbind_iobr(data1, data2, data3)
```
