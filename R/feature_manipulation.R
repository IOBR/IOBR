#' Feature Quality Control and Filtering
#'
#' @description
#' Filters features (variables) in a matrix or data frame by removing those with
#' missing values, non-numeric types, infinite values, or zero variance.
#' This is useful for preparing data for downstream statistical analyses.
#'
#' @param data A matrix or data frame containing features to filter.
#' @param feature Character vector of feature names to check. If `is_matrix = TRUE`,
#'   features are extracted from row names of the matrix.
#' @param is_matrix Logical indicating whether `data` is a gene expression matrix
#'   (features as rows, samples as columns). If `TRUE`, the matrix is transposed
#'   for processing. Default is `FALSE`.
#' @param print_result Logical indicating whether to print filtering statistics.
#'   Default is `FALSE`.
#'
#' @return Character vector of feature names that pass all quality checks.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' eset_stad <- load_data("eset_stad")
#' feas <- feature_manipulation(
#'   data = eset_stad,
#'   feature = rownames(eset_stad),
#'   is_matrix = TRUE,
#'   print_result = TRUE
#' )
#' }
feature_manipulation <- function(data,
                                 feature = NULL,
                                 is_matrix = FALSE,
                                 print_result = FALSE) {
  # Handle matrix input (features as rows)
  if (is_matrix) {
    if (!is.matrix(data) && !is.data.frame(data)) {
      cli::cli_abort("{.arg data} must be a matrix or data frame when {.code is_matrix = TRUE}")
    }
    data <- as.data.frame(t(data))
    feature <- colnames(data)
  }

  if (is.null(feature)) {
    cli::cli_abort("{.arg feature} must be specified when {.code is_matrix = FALSE}")
  }

  data <- as.data.frame(data)

  # Check which features exist in data
  valid_features <- feature[feature %in% colnames(data)]
  missing_features <- setdiff(feature, colnames(data))

  if (length(missing_features) > 0 && print_result) {
    cli::cli_alert_warning(
      "{length(missing_features)} feature{?s} not found in data: {.val {missing_features}}"
    )
  }

  if (length(valid_features) == 0) {
    cli::cli_abort("No valid features found in data")
  }

  original_count <- length(valid_features)

  # Step 1: Remove features with missing values
  if (any(is.na(data[, valid_features, drop = FALSE]))) {
    na_counts <- colSums(is.na(data[, valid_features, drop = FALSE]))
    features_with_na <- names(na_counts)[na_counts > 0]

    if (print_result) {
      cli::cli_alert_info(
        "Removing {length(features_with_na)} feature{?s} with NA values"
      )
    }

    valid_features <- setdiff(valid_features, features_with_na)
  }

  # Step 2: Remove non-numeric features
  is_numeric <- vapply(
    data[, valid_features, drop = FALSE],
    is.numeric,
    logical(1)
  )

  if (any(!is_numeric)) {
    if (print_result) {
      cli::cli_alert_info(
        "Removing {sum(!is_numeric)} non-numeric feature{?s}"
      )
    }
    valid_features <- valid_features[is_numeric]
  }

  # Step 3: Remove features with infinite values
  if (length(valid_features) > 0) {
    col_ranges <- vapply(
      data[, valid_features, drop = FALSE],
      function(x) range(x, na.rm = TRUE),
      numeric(2)
    )

    has_inf <- !is.finite(col_ranges[1, ]) | !is.finite(col_ranges[2, ])

    if (any(has_inf)) {
      if (print_result) {
        cli::cli_alert_info(
          "Removing {sum(has_inf)} feature{?s} with infinite values"
        )
      }
      valid_features <- valid_features[!has_inf]
    }
  }

  # Step 4: Remove features with zero variance
  if (length(valid_features) > 0) {
    sds <- vapply(
      data[, valid_features, drop = FALSE],
      sd,
      numeric(1),
      na.rm = TRUE
    )

    zero_sd <- sds == 0 | is.na(sds)

    if (any(zero_sd)) {
      if (print_result) {
        cli::cli_alert_info(
          "Removing {sum(zero_sd)} feature{?s} with zero variance"
        )
      }
      valid_features <- valid_features[!zero_sd]
    }
  }

  # Print summary
  if (print_result) {
    cli::cli_alert_success(
      "Retained {length(valid_features)} of {original_count} feature{?s}"
    )
  }

  valid_features
}
