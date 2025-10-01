#' Row Bind Multiple Data Sets
#'
#' Combines two or three data frames or matrices vertically using the `rbind` operation. Ensures compatibility of input data before binding.
#'
#' @param data1 Data frame or matrix. First dataset to combine.
#' @param data2 Data frame or matrix. Second dataset to combine.
#' @param data3 Data frame or matrix, optional. Additional dataset to combine. Default is NULL.
#'
#' @return Combined data frame or matrix resulting from row binding the input datasets.
#' @export
#'
#' @examples
#' data1 <- data.frame(A = 1:5, B = letters[1:5])
#' data2 <- data.frame(A = 6:10, B = letters[6:10])
#' combined_data <- rbind_iobr(data1, data2)
#' # With three datasets
#' data3 <- data.frame(A = 11:15, B = letters[11:15])
#' combined_data <- rbind_iobr(data1, data2, data3)
rbind_iobr <- function(data1, data2, data3 = NULL) {
  data2 <- assimilate_data(data_a = data1, data_b = data2)
  data.com <- rbind(data1, data2)

  if (!is.null(data3)) {
    data3 <- assimilate_data(data_a = data.com, data_b = data3)
    data.com <- rbind(data.com, data3)
  }
  return(data.com)
}
