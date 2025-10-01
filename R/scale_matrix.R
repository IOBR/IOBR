#' Scale and Manipulate a Matrix
#'
#' This function scales a given matrix, optionally applies logarithmic transformation,
#' and performs feature manipulation to remove problematic variables.
#'
#' @param matrix The matrix to be scaled and manipulated.
#' @param log2matrix Logical indicating whether to apply log2 transformation. Default is TRUE.
#' @param manipulate Logical indicating whether to perform feature manipulation. Default is TRUE.
#'
#' @return A data frame of the scaled and manipulated matrix.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("eset_gse62254", package = "IOBR")
#' eset2 <- scale_matrix(eset_gse62254, log2matrix = FALSE, manipulate = TRUE)
scale_matrix <- function(matrix, log2matrix = TRUE, manipulate = TRUE) {
  if (log2matrix) {
    matrix <- IOBR::log2eset(matrix + 1)
  }
  matrix <- as.data.frame(t(matrix))
  matrix <- scale(matrix, center = T, scale = T)
  matrix <- as.data.frame(t(matrix))

  if (manipulate) {
    feas <- feature_manipulation(data = matrix, feature = rownames(matrix), is_matrix = TRUE, print_result = TRUE)
    matrix <- matrix[rownames(matrix) %in% feas, ]
  }
  return(matrix)
}
