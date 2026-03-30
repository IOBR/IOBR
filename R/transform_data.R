#' Transform NA, Inf, or Zero Values in Data
#'
#' @description
#' Replaces NA, Inf, or zero values in specified columns of a data frame with a
#' user-defined value or the column mean.
#'
#' @param data Data frame. Input data to be transformed.
#' @param feature Character vector. Column names in `data` to apply
#'   transformation.
#' @param data_type Character. Type of value to replace: `"NA"`, `"Inf"`, or
#'   `"zero"`.
#' @param into Value to replace specified type with. Default is 0. If `"mean"`,
#'   replaces with column mean (excluding NA/Inf values as appropriate).
#'
#' @return Data frame with specified transformations applied to selected
#'   features.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' data_matrix <- data.frame(
#'   A = c(1, 2, NA, 4, Inf),
#'   B = c(Inf, 2, 3, 4, 5),
#'   C = c(0, 0, 0, 1, 2)
#' )
#'
#' # Replace NAs with 0
#' transform_data(data_matrix, feature = c("A", "B"), data_type = "NA")
#'
#' # Replace Inf values with the mean of the column
#' transform_data(data_matrix,
#'   feature = c("A", "B"),
#'   data_type = "Inf", into = "mean"
#' )
#'
#' # Replace zeros with -1 in column C
#' transform_data(data_matrix, feature = "C", data_type = "zero", into = -1)
transform_data <- function(data,
                           feature,
                           data_type = c("NA", "Inf", "zero"),
                           into = 0) {
  # Input validation
  if (!is.data.frame(data)) {
    cli::cli_abort(c(
      "Invalid input type.",
      "i" = "Expected a data frame, got {.cls {class(data)}}."
    ))
  }

  data_type <- rlang::arg_match(data_type)

  # Filter to valid features
  valid_features <- feature[feature %in% colnames(data)]
  invalid_features <- setdiff(feature, colnames(data))

  if (length(valid_features) == 0) {
    cli::cli_warn("No valid features found in data. Returning unchanged.")
    return(data)
  }

  if (length(invalid_features) > 0) {
    cli::cli_alert_warning(
      "Ignoring invalid feature{?s}: {.val {invalid_features}}"
    )
  }

  # Get column indices
  col_idx <- match(valid_features, colnames(data))

  # Define replacement function based on data_type
  replace_values <- switch(data_type,
    "NA"   = function(x) is.na(x),
    "Inf"  = function(x) is.infinite(x),
    "zero" = function(x) x == 0
  )

  # Apply transformation to each column
  for (i in seq_along(valid_features)) {
    col_name <- valid_features[i]
    j <- col_idx[i]

    mask <- replace_values(data[[j]])

    if (any(mask)) {
      if (identical(into, "mean")) {
        # Calculate mean excluding the values being replaced
        valid_vals <- data[[j]][!mask]
        replacement <- mean(valid_vals, na.rm = TRUE)
      } else {
        replacement <- into
      }
      data[mask, j] <- replacement
    }
  }

  data
}
