#' Harmonize Two Data Frames by Column Structure
#'
#' @description
#' Adds missing columns (filled with `NA`) to a secondary data frame so that
#' its column set and order match a reference data frame. This is useful when
#' combining data frames from different sources that should have the same
#' structure but may be missing some columns.
#'
#' @param data_a Data frame. Reference data frame whose column structure
#'   should be matched.
#' @param data_b Data frame. Data frame to be conformed to `data_a`.
#'
#' @return Data frame `data_b` with added missing columns (NA-filled) and
#'   reordered to match `data_a`.
#'
#' @export
#'
#' @examples
#' # Create reference data frame
#' pdata_a <- data.frame(
#'   A = 1:5, B = 2:6, C = 3:7, D = 4:8, E = 5:9
#' )
#'
#' # Create data frame with subset of columns
#' pdata_b <- data.frame(A = 1:3, C = 4:6, E = 7:9)
#'
#' # Harmonize pdata_b to match pdata_a structure
#' pdata_b_harmonized <- assimilate_data(data_a = pdata_a, data_b = pdata_b)
#' print(names(pdata_b_harmonized))  # Now has A, B, C, D, E
assimilate_data <- function(data_a, data_b) {
  # Input validation
  if (is.null(data_a) || !is.data.frame(data_a)) {
    cli::cli_abort("{.arg data_a} must be a non-null data frame.")
  }
  if (is.null(data_b) || !is.data.frame(data_b)) {
    cli::cli_abort("{.arg data_b} must be a non-null data frame.")
  }

  if (nrow(data_a) == 0) {
    cli::cli_warn("{.arg data_a} has 0 rows.")
  }
  if (nrow(data_b) == 0) {
    cli::cli_warn("{.arg data_b} has 0 rows.")
  }

  # Identify missing columns
  missing_cols <- setdiff(names(data_a), names(data_b))

  # Add missing columns with NA
  if (length(missing_cols) > 0) {
    cli::cli_alert_info(
      "Adding {length(missing_cols)} missing column{?s}: {.val {missing_cols}}"
    )

    missing_data <- stats::setNames(
      lapply(missing_cols, function(x) rep(NA, nrow(data_b))),
      missing_cols
    )
    data_b <- cbind(data_b, as.data.frame(missing_data))
  }

  # Reorder columns to match data_a
  data_b[, names(data_a), drop = FALSE]
}
