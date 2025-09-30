#' Converting mouse gene symbol to human gene symbol of expression set.
#'
#' This function converts mouse gene symbols to human gene symbols in an expression dataset. It supports using either an online resource (Ensembl) or a local dataset for conversion.
#'
#'
#' @param eset expression matrix or data frame
#' @param source default is ensembl, if an error is reported, set parameter to `local` and mus_human_gene_symbol will be use to convert gene symbols
#' @param is_matrix Boolean indicating if 'eset' is a matrix; defaults to TRUE. If FALSE, 'column_of_symbol' must be specified.
#' @param column_of_symbol default is null.Name of the column in 'eset' that contains gene symbols, if 'eset' is not a matrix. This parameter must be specified if 'is_matrix' is FALSE.
#'
#' @return Returns the expression set with human gene symbols for further analysis.
#' @export
#'
#' @examples
#' # Set the number of rows and columns
#' num_rows <- 200
#' num_cols <- 10
#'
#' # Generate sample names and mouse gene symbols
#' sample_names <- paste0("Sample", 1:num_cols)
#' mouse_genes <- paste0("MouseGene", 1:num_rows)
#'
#' # Create random data for the data frame
#' data <- matrix(runif(num_rows * num_cols), nrow = num_rows, ncol = num_cols)
#' df <- data
#'
#' rownames(df) <- anno_gc_vm32$symbol[1:200]
#' human_data <- mouse2human_eset(df, source = "local", is_matrix = TRUE)
mouse2human_eset <- function(eset, source = "ensembl", is_matrix = TRUE, column_of_symbol = NULL, verbose = FALSE) {
  if (is_matrix) {
    genes <- rownames(eset)
  } else {
    eset <- remove_duplicate_genes(eset = eset, column_of_symbol = column_of_symbol)
    genes <- rownames(eset)
  }


  if (source == "ensembl") {
    require("biomaRt")
    ensembl <- biomaRt::useEnsembl(biomart = "ensembl")
    if (verbose) print(head(listDatasets(ensembl)))
    # Basic function to convert mouse to human gene names

    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    probe_data <- getLDS(
      attributes = c("mgi_symbol"),
      filters = "mgi_symbol",
      values = genes,
      mart = mouse,
      attributesL = c("hgnc_symbol"),
      martL = human,
      uniqueRows = T
    )
    colnames(probe_data) <- c("gene_symbol_mus", "gene_symbol_human")
  } else if (source == "local") {
    probe_data <- mus_human_gene_symbol
  }

  eset <- anno_eset(eset = eset, annotation = probe_data, symbol = "gene_symbol_human", probe = "gene_symbol_mus", method = "mean")
  return(eset)
}
