# Merge Data Frames with Duplicated Column Names

Merges two data frames, resolving duplicated column names according to
user preference. Allows selection of which data frame's duplicated
columns to retain, ensuring data integrity during merging.

## Usage

``` r
merge_duplicate(
  x,
  y,
  by.x,
  by.y,
  all.x = FALSE,
  all.y = FALSE,
  all = NULL,
  choose = c("x", "y")
)
```

## Arguments

- x:

  Data frame. First data frame to merge.

- y:

  Data frame. Second data frame to merge.

- by.x:

  Character. Column name(s) in \`x\` used for merging.

- by.y:

  Character. Column name(s) in \`y\` used for merging.

- all.x:

  Logical. Include all rows from \`x\` in output. Default is \`FALSE\`.

- all.y:

  Logical. Include all rows from \`y\` in output. Default is \`FALSE\`.

- all:

  Logical or \`NULL\`. If not \`NULL\`, include all rows from both \`x\`
  and \`y\`, overriding \`all.x\` and \`all.y\`.

- choose:

  Character. Which data frame's duplicated non-joining columns to
  retain: \`"x"\` or \`"y"\`. Default is \`"x"\`.

## Value

Data frame resulting from merging \`x\` and \`y\` according to specified
parameters.

## Examples

``` r
df1 <- data.frame(ID = 1:3, Name = c("A", "B", "C"), Value = 1:3)
df2 <- data.frame(ID = 1:3, Name = c("X", "Y", "Z"), Score = 4:6)

# Merge keeping duplicated columns from x
merged_df <- merge_duplicate(df1, df2,
  by.x = "ID", by.y = "ID",
  all.x = TRUE, choose = "x"
)
#> ℹ Removing 1 duplicate column from `y`
print(merged_df)
#>   ID Name Value Score
#> 1  1    A     1     4
#> 2  2    B     2     5
#> 3  3    C     3     6

# Merge keeping duplicated columns from y
merged_df2 <- merge_duplicate(df1, df2,
  by.x = "ID", by.y = "ID",
  all = TRUE, choose = "y"
)
#> ℹ Removing 1 duplicate column from `x`
```
