#' Identify Variable Genes in Expression Data
#'
#' Identifies variable genes from a gene expression dataset using specified selection criteria. Supports multiple methods, including expression thresholding and variability estimation via median absolute deviation (MAD). Handles both count and normalized data.
#'
#' @param eset Numeric matrix. Gene expression data (genes as rows, samples as columns).
#' @param data_type Character. Type of data: "count" or "normalized". Default is "count".
#' @param methods Character vector. Methods for gene selection: "low", "mad". Default is c("low", "mad").
#' @param prop Numeric. Proportion of samples in which a gene must be expressed. Default is 0.7.
#' @param quantile Numeric vector. Quantiles for minimum MAD threshold. Default is c(0.75, 0.5, 0.25).
#' @param min.mad Numeric. Minimum allowable MAD value. Default is 0.1.
#' @param feas Character vector or NULL. Additional features to include. Default is NULL.
#'
#' @return Matrix subset of `eset` containing only variable genes identified by the specified criteria.
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' # Load expression data
#' data("eset_tme_stad", package = "IOBR")
#' # Filter variable genes
#' eset <- find_variable_genes(eset = eset_tme_stad, data_type = "normalized", methods = "mad", quantile = 0.25)
find_variable_genes <- function(eset, data_type = c("count", "normalized"), methods = c("low", "mad"), prop = 0.7,
                                quantile = c(0.75, 0.5, 0.25), min.mad = 0.1, feas = NULL) {
  eset <- as.matrix(eset)
  feas0 <- feature_manipulation(data = eset, is_matrix = TRUE)
  eset <- eset[feas0, ]

  if (data_type == "count" & c("low") %in% methods) {
    message(paste0(">>>== Genes expressed in more than", prop * 100, " percent of the samples"))
    print(table(rowSums(eset == 0) < ncol(eset) * prop))
    keep <- rowSums(eset == 0) < ncol(eset) * prop
    feas1 <- rownames(eset)[keep]
  } else {
    feas1 <- NULL
  }

  if ("mad" %in% methods) {
    eset <- log2eset(eset)
    message(paste0(">>>== min.mad = ", min.mad))
    m.mad <- apply(eset, 1, mad)
    message(paste0(">>>== Range of mad: "))
    print(range(m.mad))

    if (quantile == 0.75) index <- 4
    if (quantile == 0.50) index <- 3
    if (quantile == 0.25) index <- 2

    message(paste0(">>>== ", quantile * 100, "% of the variables will be filtered out..."))
    feas2 <- rownames(eset)[which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[index], min.mad))]
  } else {
    feas2 <- NULL
  }
  ####################################
  feas <- unique(c(feas1, feas2, feas))
  feas <- feas[!is.na(feas)]
  #####################################
  fx <- eset[rownames(eset) %in% feas, ]
  return(fx)
}
