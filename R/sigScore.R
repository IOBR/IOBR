#' Calculate Signature Score Using PCA or Mean Methods
#'
#' @description
#' Computes signature scores from gene expression data using either Principal
#' Component Analysis (PCA) or mean-based approaches.
#'
#' @param eset Normalized expression matrix with genes (signature) as rows and
#'   samples as columns.
#' @param methods Scoring method: `"PCA"` (default) for principal component 1,
#'   or `"mean"` for mean expression.
#'
#' @return Numeric vector of length `ncol(eset)`; a score summarizing the rows
#'   of `eset`.
#'
#' @author Dorothee Nickles, Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' # Load example data
#' eset_stad <- load_data("eset_stad")
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#'
#' # Get signature genes
#' signature_tme <- load_data("signature_tme")
#' genes <- signature_tme[["CD_8_T_effector"]]
#' genes <- genes[genes %in% rownames(eset)]
#'
#' # Calculate scores
#' score_pca <- sigScore(eset = eset[genes, ], methods = "PCA")
#' score_mean <- sigScore(eset = eset[genes, ], methods = "mean")
#' }
sigScore <- function(eset, methods = c("PCA", "mean")) {
  methods <- rlang::arg_match(methods)

  # Ensure numeric matrix
  eset <- as.matrix(eset)

  if (methods == "PCA") {
    # PCA-based score: PC1 weighted by correlation with mean expression
    pc <- stats::prcomp(t(eset), na.action = na.omit, scale. = TRUE)
    sigs <- pc$x[, 1] *
      sign(stats::cor(pc$x[, 1], colMeans(eset, na.rm = TRUE)))
  } else {
    # Mean-based score
    sigs <- colMeans(eset, na.rm = TRUE)
  }

  sigs
}
