


#' rbind_iobr
#'@description This function combines multiple datasets vertically using the 'rbind' operation. It takes in two required parameters, 'data1' and 'data2', and an optional parameter, 'data3'.
#'
#' @param data1 The first dataset to be combined. It should be a data.frame or a matrix.
#' @param data2 The second dataset to be combined. It should be a data.frame or a matrix.
#' @param data3  (Optional) An additional dataset to be combined. It should be a data.frame or a matrix. Default value is NULL.
#'
#' @return  The combined dataset resulting from the 'rbind' operation on the input datasets.
#' @export
#'
#' @examples
#' data1 <- data.frame(A = 1:5, B = letters[1:5])
#' data2 <- data.frame(A = 6:10, B = letters[6:10])
#' combined_data <- rbind_iobr(data1, data2)
#' 
#' # With three datasets
#' data3 <- data.frame(A = 11:15, B = letters[11:15])
#' combined_data <- rbind_iobr(data1, data2, data3)
rbind_iobr <- function(data1, data2, data3 = NULL){

  data2 <- assimilate_data(data_a = data1, data_b = data2)
  data.com <- rbind(data1, data2)

  if(!is.null(data3)){
    data3 <- assimilate_data(data_a = data.com, data_b = data3)
    data.com <- rbind(data.com, data3)
  }
  return(data.com)
}
