#' Manipulate Features in a Matrix or Data Frame
#'
#' This function performs feature manipulation on a given matrix or data frame by removing variables
#' with missing values, non-numeric variables, infinite values, and variables with zero standard deviation.
#'
#' @param data A matrix or data frame containing the features to be manipulated.
#' @param print_result Logical indicating whether to print the results of each manipulation step. Default is FALSE.
#' @param feature A vector of feature names to manipulate. Ignored if `is_matrix` is TRUE.
#' @param is_matrix Logical indicating whether the input data is in matrix format. If TRUE, feature names are extracted from column names. Default is FALSE.
#'
#' @return A vector of filtered feature names after removing problematic variables.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Load example data
#' data("eset_stad", package = "IOBR")
#' # Filter features with NA or outliers
#' feas <- feature_manipulation(data = eset_stad, feature = rownames(eset_stad), is_matrix = TRUE)
feature_manipulation <- function(data, feature, is_matrix = FALSE, print_result = FALSE) {
  if (is_matrix) {
    data <- as.data.frame(t(data))
    feature <- colnames(data)
  }

  data <- as.data.frame(data)
  # remove NA variables
  if (print_result) print(paste(">>> Is NA exist: ", sum(is.na(data[, feature]))))

  if (sum(is.na(data)) > 0) {
    nn <- as.data.frame(t(data[, feature]))
    delete_vars <- rownames(nn)[!complete.cases(nn)]
    feature <- feature[!feature %in% delete_vars]
  }
  #####################################

  # remove non-numeric variables
  if (print_result) {
    print(paste0(">>>> Is nonnumeric variables exist ? >>>>"))
    print(summary(sapply(data[, feature], mode) != "numeric"))
  }
  fea_class <- as.vector(sapply(data[, feature], mode) == "numeric")
  feature <- feature[fea_class]

  # remove infinite variables

  if (print_result) {
    print(paste0(">>>> Is -Inf variables exist ? >>>>"))
    print(summary(lapply(data[, feature], function(x) min(x)) == -Inf))
  }

  fea_class <- as.vector(lapply(data[, feature], function(x) min(x)) == -Inf)
  feature <- feature[!fea_class]

  if (print_result) {
    print(paste0(">>>> Is Inf variables exist ? >>>>"))
    print(summary(lapply(data[, feature], function(x) max(x)) == Inf))
  }
  fea_class <- as.vector(lapply(data[, feature], function(x) max(x)) == Inf)
  feature <- feature[!fea_class]

  # remove variables with same number
  sd <- apply(data[, feature], 2, function(x) sd(x) == 0)

  if (print_result) {
    print(paste0(">>> Variables with sd = 0 :  "))
    print(summary(sd))
  }

  feature <- feature[!sd]

  return(feature)
}
