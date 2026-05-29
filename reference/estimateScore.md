# estimateScore

This function reads a gene expression dataset in GCT format, calculates
enrichment scores for specific gene sets, and writes the computed scores
to an output file. It supports multiple platform types and performs
platform-specific calculations if necessary.

## Usage

``` r
estimateScore(
  input.ds,
  output.ds,
  platform = c("affymetrix", "agilent", "illumina")
)
```

## Arguments

- input.ds:

  A character string specifying the path to the input dataset file in
  GCT format. The file should have gene expression data with appropriate
  headers.

- output.ds:

  A character string specifying the path to the output dataset file,
  where the calculated scores will be written.

- platform:

  A character vector indicating the platform type. Must be one of
  "affymetrix", "agilent", or "illumina". Platform-specific calculations
  are performed based on this parameter.

## Value

This function does not return a value but writes the computed scores to
the specified output file in GCT format.

## Examples

``` r
# Simulate data
set.seed(123)
si_geneset_data <- load_data("SI_geneset")
#> ℹ Trying mirror 1/12: <https://github.com>
#> ✔ Download complete: "SI_geneset"
if (!is.null(si_geneset_data)) {
  gene_names <- unique(c(si_geneset_data[1, -1], si_geneset_data[2, -1]))
  gene_names <- gene_names[!is.na(gene_names) & gene_names != ""]
  gene_names <- head(gene_names, 500)
  n_genes <- length(gene_names)
  eset_sim <- as.data.frame(matrix(rnorm(n_genes * 3, mean = 5, sd = 1), n_genes, 3))
  rownames(eset_sim) <- gene_names
  colnames(eset_sim) <- c("Sample1", "Sample2", "Sample3")
  eset_sim <- tibble::rownames_to_column(eset_sim, var = "symbol")

  input_file <- tempfile(pattern = "estimate_", fileext = ".gct")
  output_file <- tempfile(pattern = "estimate_score_", fileext = ".gct")
  writeLines(c("#1.2", paste(nrow(eset_sim), ncol(eset_sim) - 1, sep = "\t")), input_file)
  utils::write.table(
    eset_sim,
    input_file,
    sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE
  )
  score_res <- estimateScore(
    input.ds = input_file, output.ds = output_file, platform = "affymetrix"
  )
  if (!isFALSE(score_res) && file.exists(output_file)) {
    head(read.table(output_file, skip = 2, header = TRUE, sep = "\t"))
  }
}
#> Warning: appending column names to file
#> ℹ Loading cached data: "SI_geneset"
#> 1 gene set: StromalSignature overlap=141
#> 2 gene set: ImmuneSignature overlap=141
#>            NAME   Description    Sample2    Sample3
#> 1  StromalScore  StromalScore 10.1153800 -5.4854038
#> 2   ImmuneScore   ImmuneScore 20.8432073 35.9877282
#> 3 ESTIMATEScore ESTIMATEScore 30.9585873 30.5023244
#> 4   TumorPurity   TumorPurity  0.8199163  0.8199546
```
