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
#' # Calculate scores (only if enough genes are available)
#' if (length(genes) >= 2) {
#'   score_pca <- sigScore(eset = eset[genes, ], methods = "PCA")
#'   score_mean <- sigScore(eset = eset[genes, ], methods = "mean")
#' }
#' }
sigScore <- function(eset, methods = c("PCA", "mean")) {
  methods <- rlang::arg_match(methods)

  # Ensure numeric matrix
  eset <- as.matrix(eset)

  # Check for empty matrix
  if (nrow(eset) == 0 || ncol(eset) == 0) {
    cli::cli_abort(c(
      "Expression matrix is empty.",
      "i" = "Check that input genes exist in the expression matrix.",
      "*" = "Current dimensions: {nrow(eset)} rows x {ncol(eset)} columns"
    ))
  }

  # Check for sufficient genes
  if (nrow(eset) < 2) {
    cli::cli_abort(c(
      "At least 2 genes are required for PCA method.",
      "i" = "Current number of genes: {nrow(eset)}"
    ))
  }

  # Check for samples with zero variance
  col_vars <- apply(eset, 2, stats::var, na.rm = TRUE)
  if (all(col_vars == 0, na.rm = TRUE)) {
    cli::cli_abort("All samples have zero variance.")
  }

  if (methods == "PCA") {
    # PCA-based score: PC1 weighted by correlation with mean expression
    pc <- stats::prcomp(t(eset), scale. = TRUE)
    sigs <- pc$x[, 1] *
      sign(stats::cor(pc$x[, 1], colMeans(eset, na.rm = TRUE)))
  } else {
    # Mean-based score
    sigs <- colMeans(eset, na.rm = TRUE)
  }

  sigs
}
