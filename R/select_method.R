#' Select a Signature Scoring Method Subset
#'
#' @description
#' Filters an integrated signature score matrix to retain results from a
#' specified method (PCA, ssGSEA, or zscore) and strips method suffixes from
#' column names.
#'
#' @param data Data frame or matrix. Integrated signature score matrix.
#' @param method Character. One of "PCA", "ssGSEA", or "zscore" (case-insensitive).
#'   Default is "ssGSEA".
#'
#' @return Matrix or data frame containing only the selected method's scores.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' eset_stad <- load_data("eset_stad")
#' anno_grch38 <- load_data("anno_grch38")
#' hallmark <- load_data("hallmark")
#' eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#' eset <- eset[1:5000, 1:10]
#' res <- calculate_sig_score(
#'   eset = eset,
#'   signature = hallmark[1:4],
#'   method = "integration"
#' )
#' select_method(res, method = "PCA")
select_method <- function(data, method = c("ssGSEA", "PCA", "zscore")) {
  method <- tolower(rlang::arg_match(method))

  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort("{.arg data} must be a data frame or matrix.")
  }

  patterns_to_remove <- switch(method,
    "ssgsea" = c("_PCA", "_zscore"),
    "pca"    = c("_ssGSEA", "_zscore"),
    "zscore" = c("_ssGSEA", "_PCA")
  )

  suffix_to_remove <- switch(method,
    "ssgsea" = "_ssGSEA",
    "pca"    = "_PCA",
    "zscore" = "_zscore"
  )

  for (pattern in patterns_to_remove) {
    data <- data[, !grepl(colnames(data), pattern = pattern), drop = FALSE]
  }

  colnames(data) <- gsub(colnames(data), pattern = suffix_to_remove, replacement = "")

  data
}
