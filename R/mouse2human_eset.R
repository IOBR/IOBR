#' Convert Mouse Gene Symbols to Human Gene Symbols
#'
#' @description
#' Converts mouse gene symbols to human gene symbols in an expression dataset.
#' Supports using either an online resource (Ensembl) or a local dataset for
#' conversion.
#'
#' @param eset Matrix or data frame. Expression matrix with genes in rows.
#' @param source Character. Data source for conversion: "ensembl" (online) or
#'   "local". Default is "ensembl". If Ensembl fails, use "local" which uses
#'   the internal `mus_human_gene_symbol` dataset.
#' @param is_matrix Logical. Whether `eset` is a matrix with gene symbols as
#'   row names. Default is `TRUE`. If `FALSE`, `column_of_symbol` must be
#'   specified.
#' @param column_of_symbol Character or `NULL`. Column name containing gene
#'   symbols if `eset` is not a matrix. Default is `NULL`.
#' @param verbose Logical. If `TRUE`, prints available Ensembl datasets.
#'   Default is `FALSE`.
#'
#' @return Expression set with human gene symbols.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' if (interactive()) {
#'   set.seed(123)
#'   data <- matrix(runif(50 * 3), nrow = 50, ncol = 3)
#'   rownames(data) <- c("Tpt1", "Hmgb1", "Gapdh", paste0("Gene", 4:50))
#'   colnames(data) <- paste0("Sample", 1:3)
#'   human_data <- mouse2human_eset(data, source = "local", is_matrix = TRUE)
#'   if (!is.null(human_data)) head(human_data)
#' }
mouse2human_eset <- function(eset,
                             source = c("local", "ensembl"),
                             is_matrix = TRUE,
                             column_of_symbol = NULL,
                             verbose = FALSE) {
  if (is.null(eset)) return(NULL)
  source <- rlang::arg_match(source)

  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort("{.arg eset} must be a matrix or data frame.")
  }

  if (is_matrix) {
    if (is.null(rownames(eset))) {
      cli::cli_abort("{.arg eset} must have row names when {.code is_matrix = TRUE}.")
    }
    genes <- rownames(eset)
  } else {
    if (is.null(column_of_symbol)) {
      cli::cli_abort("{.arg column_of_symbol} must be specified when {.code is_matrix = FALSE}.")
    }
    if (!column_of_symbol %in% colnames(eset)) {
      cli::cli_abort("Column {.val {column_of_symbol}} not found in {.arg eset}.")
    }
    eset <- remove_duplicate_genes(eset = eset, column_of_symbol = column_of_symbol)
    genes <- rownames(eset)
  }

  if (source == "ensembl") {
    rlang::check_installed("biomaRt")
    probe_data <- tryCatch(
      {
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl")
        if (verbose) {
          datasets <- biomaRt::listDatasets(ensembl)
          message(paste(capture.output(head(datasets)), collapse = "\n"))
        }

        human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

        probe_data <- biomaRt::getLDS(
          attributes = c("mgi_symbol"),
          filters = "mgi_symbol",
          values = genes,
          mart = mouse,
          attributesL = c("hgnc_symbol"),
          martL = human,
          uniqueRows = TRUE
        )
        colnames(probe_data) <- c("gene_symbol_mus", "gene_symbol_human")
        probe_data
      },
      error = function(e) {
        cli::cli_warn("Failed to connect to Ensembl: {e$message}")
        cli::cli_alert_info("Falling back to local database.")
        load_data("mus_human_gene_symbol")
      }
    )
  } else {
    probe_data <- load_data("mus_human_gene_symbol")
  }

  if (is.null(probe_data)) return(NULL)

  if (nrow(probe_data) == 0) {
    cli::cli_abort("No gene mappings found for the provided genes.")
  }

  eset <- anno_eset(
    eset = eset,
    annotation = probe_data,
    symbol = "gene_symbol_human",
    probe = "gene_symbol_mus",
    method = "mean"
  )

  eset
}
