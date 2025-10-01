#' Select a Signature Scoring Method Subset
#'
#' Filters an integrated signature score matrix to retain results from a specified method (PCA, ssGSEA, or zscore) and strips method suffixes from column names.
#'
#' @param data Data frame or matrix. Integrated signature score matrix.
#' @param method Character. One of "PCA", "ssGSEA", or "zscore" (case-insensitive). Default "ssGSEA".
#'
#' @return Matrix/data frame containing only the selected method's scores.
#' @export
#'
#' @examples
#' data(eset_stad, package = "IOBR")
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' eset <- scale_matrix(eset, manipulate = TRUE)
#' res <- calculate_sig_score(eset = eset, signature = signature_collection[1:4], method = "integration")
#' select_method(res, method = "PCA")
select_method <- function(data, method = "ssGSEA") {
  method <- tolower(method)

  if (method == "ssgsea") {
    data <- data[, -grep(colnames(data), pattern = "_PCA")]
    data <- data[, -grep(colnames(data), pattern = "_zscore")]
    colnames(data) <- gsub(colnames(data), pattern = "_ssGSEA", replacement = "")
  }
  if (method == "pca") {
    data <- data[, -grep(colnames(data), pattern = "_ssGSEA")]
    data <- data[, -grep(colnames(data), pattern = "_zscore")]
    colnames(data) <- gsub(colnames(data), pattern = "_PCA", replacement = "")
  }
  if (method == "zscore") {
    data <- data[, -grep(colnames(data), pattern = "_ssGSEA")]
    data <- data[, -grep(colnames(data), pattern = "_PCA")]
    colnames(data) <- gsub(colnames(data), pattern = "_zscore", replacement = "")
  }
  return(data)
}
