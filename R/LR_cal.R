#' LR_cal
#'
#' @description
#' Quantifies ligand-receptor interactions in the tumor microenvironment from bulk gene expression
#'
#' @param eset ExpressionSet object containing gene expression data.
#' @param data_type Type of input data, possible values are "count" or "tpm".
#' @param id_type Type of gene identifier, default is "ensembl".
#'
#' @references Lapuente-Santana, van Genderen, M., Hilbers, P., Finotello, F., & Eduati, F. (2021). 'Interpretable systems biomarkers predict response to immune-checkpoint inhibitors.' Patterns (New York, N.Y.), 2(8), 100293. https://doi.org/10.1016/j.patter.2021.100293
#' @return A data frame containing ligand-receptor interaction scores
#' @export
#'
#' @examples
#' data("eset_stad", package = "IOBR")
#' lr <- LR_cal(eset = eset_stad, data_type = "count", id_type = "ensembl")
LR_cal <- function(eset, data_type = c("count", "tpm"), id_type = "ensembl", cancer_type = "pancan") {
  if (!requireNamespace("easier", quietly = TRUE)) {
    stop("Package 'easier' is required but not installed. Please install it to use this function.")
  }
  # if (!requireNamespace("easier", quietly = TRUE))  BiocManager::install("easier", dependencies = FALSE)

  if (data_type == "count") {
    message(">>>=== count to TPM...")
    eset <- count2tpm(countMat = eset, idType = id_type, source = "local")
  }
  eset <- as.matrix(eset)
  # eset <- log2eset(eset)

  feas <- feature_manipulation(data = eset, feature = rownames(eset), is_matrix = TRUE)
  eset <- eset[rownames(eset) %in% feas, ]
  res <- easier::compute_LR_pairs(RNA_tpm = eset, cancer_type = cancer_type, verbose = TRUE)

  res <- rownames_to_column(res, var = "ID")
}
