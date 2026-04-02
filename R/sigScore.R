#' Calculate Signature Score Using PCA, Mean, or Z-score Methods
#'
#' @description
#' Computes signature scores from gene expression data using Principal
#' Component Analysis (PCA), mean-based, or z-score approaches.
#'
#' @param eset Normalized expression matrix with genes (signature) as rows and
#'   samples as columns.
#' @param methods Scoring method: `"PCA"` (default) for principal component 1,
#'   `"mean"` for mean expression, or `"zscore"` for z-score normalized mean.
#'
#' @return Numeric vector of length `ncol(eset)`; a score summarizing the rows
#'   of `eset`.
#'
#' @author Dorothee Nickles, Dongqiang Zeng
#' @export
#'
#' @examples
#' # Load example data
#' eset_stad <- load_data("eset_stad")
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#'
#' # Get signature genes
#' signature_tme <- load_data("signature_tme")
#' genes <- signature_tme[["CD_8_T_effector"]]
#' genes <- genes[genes %in% rownames(eset)]
#'
#' # Calculate scores (only if enough genes are available)
#' if (length(genes) >= 2) {
#'   score_pca <- sigScore(eset = eset[genes, ], methods = "PCA")
#'   score_mean <- sigScore(eset = eset[genes, ], methods = "mean")
#'   score_zscore <- sigScore(eset = eset[genes, ], methods = "zscore")
#' }
sigScore <- function(eset, methods = c("PCA", "mean", "zscore")) {
  methods <- rlang::arg_match(methods)

  eset <- as.matrix(eset)

  if (nrow(eset) == 0 || ncol(eset) == 0) {
    cli::cli_abort(c(
      "Expression matrix is empty.",
      "i" = "Check that input genes exist in the expression matrix.",
      "*" = "Current dimensions: {nrow(eset)} rows x {ncol(eset)} columns"
    ))
  }

  if (nrow(eset) < 2) {
    cli::cli_abort(c(
      "At least 2 genes are required for PCA method.",
      "i" = "Current number of genes: {nrow(eset)}"
    ))
  }

  col_vars <- apply(eset, 2, stats::var, na.rm = TRUE)
  if (all(col_vars == 0, na.rm = TRUE)) {
    cli::cli_abort("All samples have zero variance.")
  }

  if (methods == "PCA") {
    pc <- stats::prcomp(t(eset), scale. = TRUE)
    sigs <- pc$x[, 1] *
      sign(stats::cor(pc$x[, 1], colMeans(eset, na.rm = TRUE)))
  } else if (methods == "zscore") {
    eset_z <- t(scale(t(eset)))
    eset_z[is.na(eset_z)] <- 0
    sigs <- colMeans(eset_z, na.rm = TRUE)
  } else {
    sigs <- colMeans(eset, na.rm = TRUE)
  }

  sigs
}
