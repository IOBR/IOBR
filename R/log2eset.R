#' Log2 Transformation of Gene Expression Matrix
#'
#' @description
#' Determines whether a gene expression matrix requires log2 transformation
#' based on the distribution of values, and applies it if necessary. This is
#' useful for automatically detecting raw counts or linear-scale data that
#' should be log-transformed for downstream analysis.
#'
#' @param eset Numeric matrix. Gene expression data with genes as rows and
#'   samples in columns.
#'
#' @return Numeric matrix. Log2-transformed gene expression data (if
#'   transformation was needed), or the original data otherwise.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Load TCGA-STAD expression data (raw count matrix)
#' eset_stad <- load_data("eset_stad")
#'
#' # Transform count data to TPM
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#'
#' # Apply log2 transformation if needed
#' eset <- log2eset(eset)
#' }
log2eset <- function(eset) {
  # Input validation
  if (is.null(eset)) {
    cli::cli_abort("{.arg eset} cannot be NULL.")
  }

  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort(c(
      "Invalid input type.",
      "i" = "Expected a matrix or data frame, got {.cls {class(eset)}}."
    ))
  }

  # Ensure numeric matrix
  if (!is.matrix(eset)) {
    eset <- as.matrix(eset)
  }

  if (!is.numeric(eset)) {
    cli::cli_abort("Expression matrix must contain numeric values.")
  }

  # Calculate quantiles for log-judge heuristic
  qx <- stats::quantile(
    eset,
    probs = c(0, 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = TRUE
  )

  # Log-judge conditions:
  # 1. 99th percentile > 100 (large values typical of counts)
  # 2. Range > 50 AND Q1 > 0 (wide distribution, positive values)
  # 3. 0 < Q1 < 1 AND 1 < Q4 < 2 (likely already log-transformed)
  needs_log <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

  if (needs_log) {
    # Set negative values to 0 before log transformation
    eset[eset < 0] <- 0
    eset <- log2(eset + 1)
    cli::cli_alert_success("Applied log2 transformation")
  } else {
    cli::cli_alert_info(
      paste(
        "Log2 transformation not necessary",
        "(data appears to already be log-scaled)"
      )
    )
  }

  eset
}
