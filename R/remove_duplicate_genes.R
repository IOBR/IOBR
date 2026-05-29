#' Remove Duplicate Gene Symbols in Gene Expression Data
#'
#' @description
#' This function addresses duplicate gene symbols in a gene expression dataset
#' by selecting the highest-expressing instance among duplicates. Users can
#' choose between mean, standard deviation, or sum as the ranking criterion
#' for selection. This is useful for preparing data where duplicates can lead
#' to issues in downstream analyses.
#'
#' @param eset A data frame or matrix representing gene expression data, with
#'   gene symbols as one of the columns.
#' @param column_of_symbol The name of the column containing gene symbols in
#'   `eset`.
#' @param method The ranking method to use for selecting among duplicate gene
#'   symbols: `"mean"` for mean expression, `"sd"` for standard deviation,
#'   or `"sum"` for sum of expression values. Default is `"mean"`.
#'
#' @return A modified version of `eset` where duplicate gene symbols have been
#'   reduced to a single entry (the highest-ranking one). The gene symbols are
#'   set as row names in the returned data frame.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @note
#' Important: This function performs selection, not aggregation. For duplicate
#' genes, it retains only the highest-ranking instance (based on the specified
#' method) and discards others.
#'
#' @examples
#' set.seed(123)
#' test_eset <- data.frame(
#'   symbol = c("GeneA", "GeneA", "GeneB", "GeneC"),
#'   S1 = c(10, 5, 20, 15),
#'   S2 = c(12, 7, 22, 17)
#' )
#' # Remove duplicates using mean expression
#' test_eset_unique <- remove_duplicate_genes(
#'   eset = test_eset,
#'   column_of_symbol = "symbol",
#'   method = "mean"
#' )
#' print(test_eset_unique)
remove_duplicate_genes <- function(eset,
                                   column_of_symbol,
                                   method = c("mean", "sd", "sum")) {
  if (is.null(eset)) return(NULL)
  # Input validation
  if (!is.data.frame(eset) && !is.matrix(eset)) {
    cli::cli_abort(c(
      "Invalid type for {.arg eset}.",
      "i" = "Expected a data frame or matrix, got {.cls {class(eset)}}."
    ))
  }

  eset <- as.data.frame(eset)
  rownames(eset) <- NULL

  if (!column_of_symbol %in% colnames(eset)) {
    cli::cli_abort("Column {.val {column_of_symbol}} not found in eset.")
  }

  method <- rlang::arg_match(method)

  # Check for duplicates
  n_dups <- nrow(eset) - length(unique(eset[[column_of_symbol]]))

  if (n_dups == 0) {
    cli::cli_alert_info("No duplicate gene symbols found.")
    return(tibble::column_to_rownames(eset, var = column_of_symbol))
  }

  cli::cli_alert_info(
    "Found {n_dups} duplicate symbol{?s}. Using {.val {method}} for ranking."
  )

  # Get expression columns (excluding symbol column)
  expr_cols <- setdiff(colnames(eset), column_of_symbol)

  # Calculate ranking index based on method (per row)
  order_index <- switch(method,
    "mean" = rowMeans(eset[, expr_cols, drop = FALSE], na.rm = TRUE),
    "sd" = apply(eset[, expr_cols, drop = FALSE], 1, stats::sd, na.rm = TRUE),
    "sum" = rowSums(eset[, expr_cols, drop = FALSE], na.rm = TRUE)
  )

  # Sort by ranking (descending) and keep first occurrence
  eset <- eset[order(order_index, decreasing = TRUE), , drop = FALSE]
  eset <- dplyr::distinct(eset, !!rlang::sym(column_of_symbol),
    .keep_all = TRUE
  )

  cli::cli_alert_success(
    "Reduced to {nrow(eset)} unique gene{?s}"
  )

  tibble::column_to_rownames(eset, var = column_of_symbol)
}
