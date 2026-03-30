#' Batch Correlation Analysis
#'
#' @description
#' Performs correlation analysis between a target variable and multiple feature
#' variables. Computes correlation coefficients, p-values, and adjusts for
#' multiple testing using the Benjamini-Hochberg method.
#'
#' @param data Data frame containing the target and feature variables.
#' @param target Character string specifying the name of the target variable.
#' @param feature Character vector specifying the names of feature variables to
#'   correlate with the target.
#' @param method Character string specifying the correlation method. Options are
#'   `"spearman"`, `"pearson"`, or `"kendall"`. Default is `"spearman"`.
#'
#' @return Tibble containing the following columns for each feature:
#' \describe{
#'   \item{sig_names}{Feature name}
#'   \item{p.value}{Raw p-value}
#'   \item{statistic}{Correlation coefficient}
#'   \item{p.adj}{Adjusted p-value (Benjamini-Hochberg)}
#'   \item{log10pvalue}{Negative log10-transformed p-value}
#'   \item{stars}{Significance stars: **** p<0.0001, *** p<0.001,
#'     ** p<0.01, * p<0.05, + p<0.5}
#' }
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' # Load TCGA-STAD signature data
#' sig_stad <- load_data("sig_stad")
#'
#' # Perform batch correlation
#' results <- batch_cor(
#'   data = sig_stad,
#'   target = "CD_8_T_effector",
#'   feature = colnames(sig_stad)[69:ncol(sig_stad)]
#' )
#' head(results)
#' }
batch_cor <- function(data,
                      target,
                      feature,
                      method = c("spearman", "pearson", "kendall")) {
  # Input validation
  if (is.null(data) || !is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a non-null data frame.")
  }
  if (nrow(data) == 0) {
    cli::cli_abort("{.arg data} has no rows.")
  }
  if (!target %in% colnames(data)) {
    cli::cli_abort("Target {.val {target}} not found in data columns.")
  }
  if (!is.character(feature) || length(feature) == 0) {
    cli::cli_abort("{.arg feature} must be a non-empty character vector.")
  }

  method <- rlang::arg_match(method)

  # Filter valid features
  valid_features <- feature[feature %in% colnames(data)]
  invalid_features <- setdiff(feature, colnames(data))

  if (length(invalid_features) > 0) {
    cli::cli_alert_warning(paste(
      "Ignoring {length(invalid_features)} invalid feature{?s}:",
      "{.val {invalid_features}}"
    ))
  }

  if (length(valid_features) == 0) {
    cli::cli_abort("No valid features found in data.")
  }

  valid_features <- setdiff(valid_features, target)

  # Remove rows with NA in target
  data <- data[!is.na(data[[target]]), , drop = FALSE]

  if (nrow(data) < 3) {
    cli::cli_abort("Insufficient data after removing NA values (need >= 3).")
  }

  # Filter features with zero variance
  valid_features <- valid_features[
    vapply(
      data[, valid_features, drop = FALSE],
      function(x) stats::sd(x, na.rm = TRUE) > 0,
      logical(1)
    )
  ]

  if (length(valid_features) == 0) {
    cli::cli_abort("All features have zero variance.")
  }

  cli::cli_alert_info(
    "Computing {method} correlation for {length(valid_features)} feature{?s}"
  )

  # Vectorized correlation analysis
  results <- vapply(valid_features, function(feat) {
    test <- stats::cor.test(data[[feat]], data[[target]],
      method = method, use = "complete.obs"
    )
    p_val <- exact_pvalue(data[[feat]], data[[target]], method)
    c(p.value = p_val, estimate = test$estimate)
  }, numeric(2))

  # Build results data frame
  cc <- data.frame(
    sig_names = valid_features,
    p.value = results["p.value", ],
    statistic = results["estimate", ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  cc$p.adj <- stats::p.adjust(cc$p.value, method = "BH")
  cc$log10pvalue <- -log10(cc$p.value)
  cc$stars <- cut(cc$p.value,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
    labels = c("****", "***", "**", "*", "+", "")
  )

  cc <- cc[order(cc$p.value), , drop = FALSE]
  rownames(cc) <- NULL

  cli::cli_alert_success("Correlation analysis complete")

  tibble::as_tibble(cc)
}


#' Calculate Exact P-Value for Correlation
#'
#' @description
#' Computes the exact p-value for the correlation between two numeric variables
#' using a specified correlation method.
#'
#' @param x Numeric vector representing the first variable.
#' @param y Numeric vector representing the second variable.
#' @param method Character string specifying the correlation method:
#'   `"spearman"`, `"pearson"`, or `"kendall"`.
#'
#' @return Numeric value representing the exact p-value.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' sig_stad <- load_data("sig_stad")
#' p_val <- exact_pvalue(
#'   x = sig_stad$CD8.T.cells,
#'   y = sig_stad$CD_8_T_effector,
#'   method = "spearman"
#' )
#' print(p_val)
#' }
exact_pvalue <- function(x, y, method) {
  l <- sum(!is.na(x) & !is.na(y))

  if (l < 3) {
    return(NA_real_)
  }

  r <- stats::cor(x = x, y = y, method = method, use = "complete.obs")

  # t-statistic for correlation
  t_stat <- r * sqrt((l - 2) / (1 - r^2))

  # Two-tailed p-value
  p <- 2 * stats::pt(q = abs(t_stat), df = l - 2, lower.tail = FALSE)

  p
}
