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
#' \donttest{
#' # Load and annotate expression data
#' eset_stad <- load_data("eset_stad")
#' anno_rnaseq <- load_data("anno_rnaseq")
#' eset_stad <- anno_eset(eset = eset_stad, annotation = anno_rnaseq)
#' eset_stad <- tibble::rownames_to_column(as.data.frame(eset_stad), var = "id")
#'
#' # Create duplicate gene names for demonstration
#' eset_stad[2:3, "id"] <- "MT-CO1"
#'
#' # Check duplicates before
#' sum(duplicated(eset_stad$id))
#'
#' # Remove duplicates using mean expression as ranking criterion
#' eset_stad <- remove_duplicate_genes(
#'   eset = eset_stad,
#'   column_of_symbol = "id",
#'   method = "mean"
#' )
#'
#' # Check duplicates after
#' sum(duplicated(rownames(eset_stad)))
#' }
remove_duplicate_genes <- function(eset,
                                    column_of_symbol,
                                    method = c("mean", "sd", "sum")) {

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

  # Calculate ranking index based on method
  order_index <- switch(method,
    "mean" = vapply(eset[, expr_cols, drop = FALSE], mean, numeric(1), na.rm = TRUE),
    "sd"   = vapply(eset[, expr_cols, drop = FALSE], stats::sd, numeric(1), na.rm = TRUE),
    "sum"  = rowSums(eset[, expr_cols, drop = FALSE], na.rm = TRUE)
  )

  # Sort by ranking (descending) and keep first occurrence
  eset <- eset[order(order_index, decreasing = TRUE), , drop = FALSE]
  eset <- dplyr::distinct(eset, !!rlang::sym(column_of_symbol), .keep_all = TRUE)

  cli::cli_alert_success(
    "Reduced to {nrow(eset)} unique gene{?s}"
  )

  tibble::column_to_rownames(eset, var = column_of_symbol)
}
