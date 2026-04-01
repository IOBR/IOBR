#' Row Bind Multiple Data Sets
#'
#' @description
#' Combines two or three data frames or matrices vertically using `rbind`.
#' Ensures compatibility of input data before binding by aligning columns.
#'
#' @param data1 Data frame or matrix. First dataset to combine.
#' @param data2 Data frame or matrix. Second dataset to combine.
#' @param data3 Data frame or matrix or `NULL`. Optional third dataset.
#'   Default is `NULL`.
#'
#' @return Combined data frame resulting from row binding the input datasets.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data1 <- data.frame(A = 1:5, B = letters[1:5])
#' data2 <- data.frame(A = 6:10, B = letters[6:10])
#' combined_data <- rbind_iobr(data1, data2)
#'
#' # With three datasets
#' data3 <- data.frame(A = 11:15, B = letters[11:15])
#' combined_data <- rbind_iobr(data1, data2, data3)
rbind_iobr <- function(data1, data2, data3 = NULL) {
  if (is.null(data1) || is.null(data2)) {
    cli::cli_abort("{.arg data1} and {.arg data2} must not be NULL.")
  }

  if (!is.data.frame(data1) && !is.matrix(data1)) {
    cli::cli_abort("{.arg data1} must be a data frame or matrix.")
  }
  if (!is.data.frame(data2) && !is.matrix(data2)) {
    cli::cli_abort("{.arg data2} must be a data frame or matrix.")
  }

  data2 <- assimilate_data(data_a = data1, data_b = data2)
  data.com <- rbind(data1, data2)

  if (!is.null(data3)) {
    if (!is.data.frame(data3) && !is.matrix(data3)) {
      cli::cli_abort("{.arg data3} must be a data frame, matrix, or NULL.")
    }
    data3 <- assimilate_data(data_a = data.com, data_b = data3)
    data.com <- rbind(data.com, data3)
  }

  data.com
}
