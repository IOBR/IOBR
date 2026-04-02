#' Scale and Manipulate a Matrix
#'
#' @description
#' Scales a gene expression matrix, optionally applies logarithmic
#' transformation, and performs feature manipulation to remove problematic
#' variables (e.g., genes with zero variance or missing values).
#'
#' @param matrix Numeric matrix with genes as rows and samples as columns.
#' @param log2matrix Logical indicating whether to apply log2 transformation
#'   using [log2eset()]. Default is `TRUE`.
#' @param manipulate Logical indicating whether to perform feature manipulation
#'   to remove problematic features. Default is `TRUE`.
#'
#' @return A scaled matrix (genes as rows, samples as columns).
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' eset_gse62254 <- load_data("eset_gse62254")
#' eset2 <- scale_matrix(eset_gse62254, log2matrix = FALSE, manipulate = TRUE)
scale_matrix <- function(matrix, log2matrix = TRUE, manipulate = TRUE) {
  # Input validation
  if (!is.matrix(matrix) && !is.data.frame(matrix)) {
    cli::cli_abort(c(
      "Invalid input type for {.arg matrix}.",
      "i" = "Expected a matrix or data frame, got {.cls {class(matrix)}}."
    ))
  }

  # Ensure numeric matrix
  if (!is.matrix(matrix)) {
    matrix <- as.matrix(matrix)
  }

  if (!is.numeric(matrix)) {
    cli::cli_abort("Matrix must contain numeric values.")
  }

  # Apply log2 transformation if requested
  if (log2matrix) {
    matrix <- log2eset(matrix)
  }

  # Scale: center and scale by column (sample), then transpose back
  # Transpose so samples become rows for scale()
  matrix <- as.data.frame(t(matrix))
  matrix <- scale(matrix, center = TRUE, scale = TRUE)
  matrix <- as.data.frame(t(matrix))

  # Feature manipulation to remove problematic features
  if (manipulate) {
    feas <- feature_manipulation(
      data = matrix,
      feature = rownames(matrix),
      is_matrix = TRUE,
      print_result = TRUE
    )
    matrix <- matrix[rownames(matrix) %in% feas, , drop = FALSE]
    cli::cli_alert_info("Retained {length(feas)} features after manipulation")
  }

  as.matrix(matrix)
}
