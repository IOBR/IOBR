# Save Signature Data to File

Saves signature data to a specified file format, supporting CSV or
RData. Handles single signatures or lists of signatures, converting them
to a data frame for storage.

## Usage

``` r
output_sig(signatures, format = c("csv", "rdata"), file.name)
```

## Arguments

- signatures:

  Signature data: a list or single string of signatures.

- format:

  Character. Output format: "csv" or "rdata". Default is "csv".

- file.name:

  Character. Name of the output file without extension.

## Value

Data frame containing the processed signature data, also saved to the
specified file.

## Author

Dongqiang Zeng

## Examples

``` r
# Simulate signatures
set.seed(123)
sim_sigs <- list(
  Signature1 = paste0("Gene", 1:50),
  Signature2 = paste0("Gene", 51:100),
  Signature3 = paste0("Gene", 101:150)
)
tmpfile <- tempfile(fileext = ".csv")
output_sig(
  signatures = sim_sigs, format = "csv",
  file.name = tools::file_path_sans_ext(tmpfile)
)
#> ✔ Signature data saved to /tmp/RtmpRGChSj/file1ccb7834e87b.csv
#>    Signature1 Signature2 Signature3
#> 1       Gene1     Gene51    Gene101
#> 2       Gene2     Gene52    Gene102
#> 3       Gene3     Gene53    Gene103
#> 4       Gene4     Gene54    Gene104
#> 5       Gene5     Gene55    Gene105
#> 6       Gene6     Gene56    Gene106
#> 7       Gene7     Gene57    Gene107
#> 8       Gene8     Gene58    Gene108
#> 9       Gene9     Gene59    Gene109
#> 10     Gene10     Gene60    Gene110
#> 11     Gene11     Gene61    Gene111
#> 12     Gene12     Gene62    Gene112
#> 13     Gene13     Gene63    Gene113
#> 14     Gene14     Gene64    Gene114
#> 15     Gene15     Gene65    Gene115
#> 16     Gene16     Gene66    Gene116
#> 17     Gene17     Gene67    Gene117
#> 18     Gene18     Gene68    Gene118
#> 19     Gene19     Gene69    Gene119
#> 20     Gene20     Gene70    Gene120
#> 21     Gene21     Gene71    Gene121
#> 22     Gene22     Gene72    Gene122
#> 23     Gene23     Gene73    Gene123
#> 24     Gene24     Gene74    Gene124
#> 25     Gene25     Gene75    Gene125
#> 26     Gene26     Gene76    Gene126
#> 27     Gene27     Gene77    Gene127
#> 28     Gene28     Gene78    Gene128
#> 29     Gene29     Gene79    Gene129
#> 30     Gene30     Gene80    Gene130
#> 31     Gene31     Gene81    Gene131
#> 32     Gene32     Gene82    Gene132
#> 33     Gene33     Gene83    Gene133
#> 34     Gene34     Gene84    Gene134
#> 35     Gene35     Gene85    Gene135
#> 36     Gene36     Gene86    Gene136
#> 37     Gene37     Gene87    Gene137
#> 38     Gene38     Gene88    Gene138
#> 39     Gene39     Gene89    Gene139
#> 40     Gene40     Gene90    Gene140
#> 41     Gene41     Gene91    Gene141
#> 42     Gene42     Gene92    Gene142
#> 43     Gene43     Gene93    Gene143
#> 44     Gene44     Gene94    Gene144
#> 45     Gene45     Gene95    Gene145
#> 46     Gene46     Gene96    Gene146
#> 47     Gene47     Gene97    Gene147
#> 48     Gene48     Gene98    Gene148
#> 49     Gene49     Gene99    Gene149
#> 50     Gene50    Gene100    Gene150
```
