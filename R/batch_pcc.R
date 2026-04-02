#' Batch Calculation of Partial Correlation Coefficients
#'
#' @description
#' Computes partial correlation coefficients between multiple features and a
#' target variable while controlling for an interference (confounding) variable.
#' Adjusts p-values for multiple testing using the Benjamini-Hochberg method.
#'
#' @param input Data frame containing feature variables, target variable, and
#'   interference variable.
#' @param interferenceid Character string specifying the column name of the
#'   interference (confounding) variable to control for.
#' @param target Character string specifying the column name of the target
#'   variable.
#' @param features Character vector specifying the column names of feature
#'   variables to correlate with the target.
#' @param method Character string specifying the correlation method. Options are
#'   `"pearson"`, `"spearman"`, or `"kendall"`. Default is `"pearson"`.
#'
#' @return Tibble containing the following columns for each feature:
#' \describe{
#'   \item{sig_names}{Feature name}
#'   \item{p.value}{Raw p-value}
#'   \item{statistic}{Partial correlation coefficient}
#'   \item{p.adj}{Adjusted p-value (Benjamini-Hochberg method)}
#'   \item{log10pvalue}{Negative log10-transformed p-value}
#'   \item{stars}{Significance stars: **** p.adj<0.0001, *** p.adj<0.001,
#'     ** p.adj<0.01, * p.adj<0.05, + p.adj<0.5}
#' }
#'
#' @author Rongfang Shen
#' @export
#'
#' @examples
#' # Load TCGA-STAD signature data
#' sig_stad <- load_data("sig_stad")
#'
#' # Calculate partial correlations controlling for tumor purity
#' res <- batch_pcc(
#'   input = sig_stad,
#'   interferenceid = "TumorPurity_estimate",
#'   target = "Pan_F_TBRs",
#'   method = "pearson",
#'   features = colnames(sig_stad)[70:ncol(sig_stad)]
#' )
#' head(res)
batch_pcc <- function(input,
                      interferenceid,
                      target,
                      features,
                      method = c("pearson", "spearman", "kendall")) {
  # Input validation
  if (is.null(input) || !is.data.frame(input)) {
    cli::cli_abort("{.arg input} must be a non-null data frame.")
  }
  if (nrow(input) == 0) {
    cli::cli_abort("{.arg input} has no rows.")
  }
  if (!interferenceid %in% colnames(input)) {
    cli::cli_abort(
      "Interference variable {.val {interferenceid}} not found in data."
    )
  }
  if (!target %in% colnames(input)) {
    cli::cli_abort("Target variable {.val {target}} not found in data.")
  }
  if (!is.character(features) || length(features) == 0) {
    cli::cli_abort("{.arg features} must be a non-empty character vector.")
  }

  method <- rlang::arg_match(method)

  # Filter valid features
  features <- setdiff(features, c(interferenceid, target))
  valid_features <- features[features %in% colnames(input)]
  invalid_features <- setdiff(features, colnames(input))

  if (length(invalid_features) > 0) {
    cli::cli_alert_warning(
      "Ignoring {length(invalid_features)} invalid feature{?s}:",
      " {.val {invalid_features}}}"
    )
  }

  if (length(valid_features) == 0) {
    cli::cli_abort("No valid features found in data.")
  }

  cli::cli_alert_info(
    "Computing {method} partial correlation for",
    " {length(valid_features)} feature{?s}"
  )

  # Partial correlation test function
  pcor_test <- function(x, y, z, method) {
    dat <- stats::na.omit(cbind(x, y, z))
    if (nrow(dat) < 4) {
      return(c(estimate = NA_real_, p.value = NA_real_))
    }

    r_mat <- stats::cor(dat, method = method)
    rxy <- r_mat[1, 2]
    rxz <- r_mat[1, 3]
    ryz <- r_mat[2, 3]

    # Partial correlation formula
    rho <- (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))

    n <- nrow(dat)
    if (abs(rho) >= 1) {
      return(c(estimate = rho, p.value = 0))
    }

    tstat <- rho * sqrt((n - 3) / (1 - rho^2))
    pval <- 2 * stats::pt(abs(tstat), df = n - 3, lower.tail = FALSE)

    c(estimate = rho, p.value = pval)
  }

  # Vectorized partial correlation calculation
  results <- vapply(valid_features, function(feat) {
    pcor_test(
      input[[feat]],
      input[[target]],
      input[[interferenceid]],
      method = method
    )
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
  cc$stars <- cut(cc$p.adj,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
    labels = c("****", "***", "**", "*", "+", "")
  )

  cc <- cc[order(cc$p.value), , drop = FALSE]

  cli::cli_alert_success("Partial correlation analysis complete")

  tibble::as_tibble(cc)
}
