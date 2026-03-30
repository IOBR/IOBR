#' Check Integrity and Outliers of Expression Set
#'
#' @description
#' Performs quality checks on an expression matrix to identify missing values,
#' infinite values, and features with zero variance. Issues warnings when
#' potential problems are detected that may affect downstream analyses.
#'
#' @param eset Expression matrix or data frame with genes/features in rows and
#'   samples in columns.
#' @param print_result Logical indicating whether to print detailed check
#'   results to the console. Default is `FALSE`.
#' @param estimate_sd Logical indicating whether to check for features with
#'   zero standard deviation. Default is `FALSE`.
#'
#' @return Invisibly returns `NULL`. Function is called for its side effects
#'   (printing messages and issuing warnings).
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' # Load TCGA-STAD expression data
#' eset_stad <- load_data("eset_stad")
#'
#' # Convert counts to TPM
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#'
#' # Check expression set integrity
#' check_eset(eset)
#'
#' # Check with detailed output
#' check_eset(eset, print_result = TRUE, estimate_sd = TRUE)
check_eset <- function(eset, print_result = FALSE, estimate_sd = FALSE) {
  # Input validation
  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort(c(
      "Invalid input type for {.arg eset}.",
      "i" = "Expected a matrix or data frame, got {.cls {class(eset)}}."
    ))
  }

  # Check for NA values
  na_count <- sum(is.na(eset))

  if (print_result) {
    cli::cli_alert_info("Checking for NA values: {na_count} found")
  }

  if (na_count > 0) {
    cli::cli_warn(c(
      paste("There are", na_count, "missing values in the matrix."),
      "i" = "This may affect score calculation.",
      "*" = paste(
        "Set adjust_eset = TRUE to handle missing values",
        "automatically."
      ),
    ))
  }

  # Check for infinite values (vectorized)
  eset_range <- range(eset, na.rm = TRUE)
  has_neg_inf <- !is.finite(eset_range[1]) && eset_range[1] < 0
  has_pos_inf <- !is.finite(eset_range[2]) && eset_range[2] > 0

  if (print_result) {
    cli::cli_alert_info(paste(
      "Checking for -Inf values:",
      sum(is.infinite(eset) & eset < 0, na.rm = TRUE),
      "found"
    ))
    cli::cli_alert_info(paste(
      "Checking for +Inf values:",
      sum(is.infinite(eset) & eset > 0, na.rm = TRUE),
      "found"
    ))
  }

  if (has_neg_inf || has_pos_inf) {
    cli::cli_warn(c(
      "Infinite values detected in the matrix.",
      "i" = "This may affect score calculation.",
      "*" = paste(
        "Set adjust_eset = TRUE to handle infinite values",
        "automatically."
      )
    ))
  }

  # Check for zero variance features
  if (estimate_sd) {
    row_sds <- apply(eset, 1, stats::sd, na.rm = TRUE)
    zero_sd_count <- sum(row_sds == 0, na.rm = TRUE)

    if (print_result) {
      cli::cli_alert_info("Features with zero variance: {zero_sd_count}")
    }

    if (zero_sd_count > 0) {
      cli::cli_warn(c(
        "{zero_sd_count} feature{?s} {?has/have} zero variance.",
        "i" = "Zero-variance features may affect score calculation.",
        "*" = "Set {.code adjust_eset = TRUE} to remove them automatically."
      ))
    }
  }

  invisible(NULL)
}
