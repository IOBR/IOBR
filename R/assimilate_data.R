#' Harmonize Two Data Frames by Column Structure
#'
#' Adds missing columns (filled with NA) to a secondary data frame so that its column set and order match a reference data frame.
#'
#' @param data_a Data frame. Reference data frame whose column structure should be matched.
#' @param data_b Data frame. Data frame to be conformed to `data_a`.
#'
#' @return Data frame `data_b` with added missing columns (NA-filled) and reordered to match `data_a`.
#' @export
#'
#' @examples
#' pdata_a <- as.data.frame(array(3, c(5, 9)))
#' colnames(pdata_a) <- LETTERS[1:dim(pdata_a)[2]]
#' pdata_b <- as.data.frame(array(2, c(7, 4)))
#' colnames(pdata_b) <- c("A", "C", "E", "F")
#' assimilate_data(data_a = pdata_a, data_b = pdata_b)
#'
assimilate_data <- function(data_a, data_b) {
  data_a <- as.data.frame(data_a)
  data_b <- as.data.frame(data_b)
  miss <- as.character(colnames(data_a)[!colnames(data_a) %in% colnames(data_b)])
  ncolumns <- length(miss)
  nrows <- dim(data_b)[1]
  missdata <- as.data.frame(array(NA, c(nrows, length(miss))))
  colnames(missdata) <- miss
  data_b <- as.data.frame(cbind(data_b, missdata))
  data_b <- data_b[, match(colnames(data_a), colnames(data_b))]
  return(data_b)
}
