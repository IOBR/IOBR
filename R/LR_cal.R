





#' LR_cal
#'
#' @description
#' Quantifies ligand-receptor interactions in the tumor microenvironment from bulk gene expression
#'
#' @param eset ExpressionSet object containing gene expression data.
#' @param data_type Type of input data, possible values are "count" or "tpm".
#' @param id_type Type of gene identifier, default is "ensembl".
#'
#' @return
#' @export
#'
#' @examples
#' data("eset_stad", package = "IOBR")
#' lr <- LR_cal(eset = eset_stad, data_type = "count", id_type = "ensembl")
#'
LR_cal <- function(eset, data_type = c("count", "tpm"), id_type = "ensembl"){

  if (!requireNamespace("easier", quietly = TRUE))  BiocManager::install("easier", dependencies = FALSE)

  if(data_type=="count"){

    message(">>>=== count to TPM...")
    eset <- count2tpm(countMat = eset, idType = id_type, source = "local")
  }
  eset <- as.matrix(eset)
  # eset <- log2eset(eset)

  feas <- feature_manipulation(data = eset, feature = rownames(eset), is_matrix = TRUE)
  eset <- eset[rownames(eset)%in%feas, ]
  res <- easier::compute_LR_pairs(RNA_tpm = eset, cancer_type = "pancan", verbose = TRUE)

  res <- rownames_to_column(res, var = "ID")

}




