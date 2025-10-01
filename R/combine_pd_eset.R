#' Combine Phenotype Data and Expression Set
#'
#' @description
#' Merges phenotype data with an expression matrix by matching sample IDs. Optionally
#' filters features, applies feature manipulation, and scales expression data before
#' combining.
#'
#' @param eset Expression matrix with genes/features in rows and samples in columns.
#' @param pdata Data frame containing phenotype/clinical data.
#' @param id_pdata Character string specifying the column name in \code{pdata} containing
#'   sample identifiers. Default is \code{"ID"}.
#' @param feas Character vector specifying features to include from \code{eset}. If
#'   \code{NULL}, all features are used. Default is \code{NULL}.
#' @param feature_manipulation Logical indicating whether to apply feature manipulation
#'   to filter valid features. Default is \code{TRUE}.
#' @param scale Logical indicating whether to scale (standardize) expression data.
#'   Default is \code{TRUE}.
#' @param choose_who_when_duplicate Character string specifying which data to prefer when
#'   duplicate columns exist. Options are \code{"eset"} or \code{"pdata"}. Default is
#'   \code{"eset"}.
#'
#' @return Data frame combining phenotype data and (transposed) expression data, with
#'   samples in rows and features/phenotypes in columns.
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD expression data
#' data("eset_stad", package = "IOBR")
#' # Convert to TPM
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' colnames(eset) <- substring(colnames(eset), 1, 12)
#' # Load phenotype data
#' data("sig_stad", package = "IOBR")
#' # Combine phenotype and expression data
#' input <- combine_pd_eset(eset = eset, pdata = sig_stad,
#'                          feas = unique(unlist(signature_collection)))
combine_pd_eset <- function(eset, pdata, id_pdata = "ID", feas = NULL, feature_manipulation = TRUE, scale = TRUE, choose_who_when_duplicate = "eset") {
  if (!is.null(feas)) {
    eset <- eset[rownames(eset) %in% feas, ]
  }

  feas_filter <- feature_manipulation(data = eset, feature = rownames(eset), is_matrix = T)

  eset <- eset[rownames(eset) %in% feas_filter, ]
  eset <- t(eset)
  if (scale) eset <- scale(eset, scale = T, center = T)

  eset <- rownames_to_column(as.data.frame(eset), var = "ID")

  colnames(pdata)[which(colnames(pdata) == id_pdata)] <- "ID"

  if (choose_who_when_duplicate == "eset") {
    choose <- "x"
  } else {
    choose <- "y"
  }
  pd_eset <- merge_duplicate(pdata, eset, by.x = "ID", by.y = "ID", choose = choose, all = FALSE)
  return(pd_eset)
}
