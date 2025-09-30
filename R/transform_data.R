#' Transform features with Inf or NA to other number
#'
#' @param data A data frame containing the data to be transformed.
#' @param feature A vector of column names in 'data' where transformations will be applied.
#' @param data_type The type of data to be transformed: "NA", "Inf", or "zero".
#' @param into The value to replace 'data_type' with; default is 0. If 'into' is "mean", it replaces with the column mean.
#'
#' @return Returns the data frame with specified transformations applied to the indicated features.
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' data_matrix <- data.frame(A = c(1, 2, NA, 4, Inf), B = c(Inf, 2, 3, 4, 5), C = c(0, 0, 0, 1, 2))
#' # Replace NAs with 0
#' transform_data(data_matrix, feature = c("A", "B"), data_type = "NA")
#' # Replace Infinite values with the mean of the column
#' transform_data(data_matrix, feature = c("A", "B"), data_type = "Inf", into = "mean")
#' # Replace zeros with -1 in column C
#' transform_data(data_matrix, feature = "C", data_type = "zero", into = -1)
transform_data <- function(data, feature, data_type = c("NA", "Inf", "zero"), into = 0) {
  feature <- colnames(data)[colnames(data) %in% feature]

  if (data_type == "NA") {
    for (i in 1:length(feature)) {
      var <- feature[i]
      j <- which(colnames(data) == var)

      if (into == "mean") {
        data[is.na(data[, var]), j] <- mean(data[, j], na.rm = TRUE)
      } else {
        data[is.na(data[, var]), j] <- into
      }
    }
  } else if (data_type == "Inf") {
    for (i in 1:length(feature)) {
      var <- feature[i]
      j <- which(colnames(data) == var)

      if (into == "mean") {
        data[is.infinite(data[, var]), j] <- mean(data[, j], na.rm = TRUE)
      } else {
        data[is.infinite(data[, var]), j] <- into
      }
    }
  } else if (data_type == "zero") {
    for (i in 1:length(feature)) {
      var <- feature[i]
      j <- which(colnames(data) == var)
      data[data[, var] == 0, j] <- into
    }
  }
  return(data)
}
