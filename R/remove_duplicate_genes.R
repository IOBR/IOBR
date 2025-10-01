#' Remove Duplicate Gene Symbols in Gene Expression Data
#'
#' This function addresses duplicate gene symbols in a gene expression dataset by aggregating
#' the expression data for duplicate entries. Users can choose between mean, standard deviation, or sum
#' for aggregation. This is useful for preparing data where duplicates can lead to issues in downstream analyses.
#'
#' @param eset A data frame or matrix representing gene expression data, with gene symbols as one of the columns.
#' @param column_of_symbol The name of the column containing gene symbols in `eset`.
#' @param method The aggregation method to apply for duplicate gene symbols: "mean" for averaging,
#'        "sd" for standard deviation, or "sum" for the sum of values. Default is "mean".
#'
#' @return A modified version of `eset` where duplicate gene symbols have been aggregated according to the specified method.
#'         The gene symbols are set as row names in the returned data frame or matrix.
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#'
#' # loading eset
#' data("eset_stad", package = "IOBR")
#' # annotation
#' eset_stad <- anno_eset(eset = eset_stad, annotation = anno_rnaseq)
#' eset_stad <- rownames_to_column(eset_stad, var = "id")
#'
#' # Creating duplicate gene names
#' eset_stad[2:3, "id"] <- "MT-CO1"
#' # Counting the number of identical names
#' summary(duplicated(eset_stad$id))
#' # De-duplication of rows with the same gene name using the average value
#' eset_stad <- remove_duplicate_genes(eset = eset_stad, column_of_symbol = "id", method = "mean")
#' summary(duplicated(eset_stad$id))
remove_duplicate_genes <- function(eset, column_of_symbol, method = "mean") {
  eset <- as.data.frame(eset)
  rownames(eset) <- NULL

  dups <- dim(eset)[1] - length(unique(eset[, column_of_symbol]))

  if (dups == 0) {
    eset <- tibble::column_to_rownames(eset, var = column_of_symbol)
    return(eset)
  } else {
    if (method == "mean") {
      order_index <- apply(eset[, setdiff(colnames(eset), column_of_symbol)], 1, function(x) mean(x, na.rm = T))
      eset <- eset[order(order_index, decreasing = T), ]
      eset <- eset %>%
        dplyr::distinct(!!sym(column_of_symbol), .keep_all = TRUE) %>%
        tibble::column_to_rownames(., var = column_of_symbol)
      return(eset)
    } else if (method == "sd") {
      order_index <- apply(eset[, setdiff(colnames(eset), column_of_symbol)], 1, function(x) sd(x, na.rm = T))
      eset <- eset[order(order_index, decreasing = T), ]
      eset <- eset %>%
        distinct(!!sym(column_of_symbol), .keep_all = TRUE) %>%
        tibble::column_to_rownames(., var = column_of_symbol)
      return(eset)
    } else if (method == "sum") {
      order_index <- apply(eset[, setdiff(colnames(eset), column_of_symbol)], 1, function(x) sum(x, na.rm = T))
      eset <- eset[order(order_index, decreasing = T), ]
      eset <- eset %>%
        distinct(!!sym(column_of_symbol), .keep_all = TRUE) %>%
        tibble::column_to_rownames(., var = column_of_symbol)
      return(eset)
    }
  }
}
