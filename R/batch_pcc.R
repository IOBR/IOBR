#' Batch Calculation of Partial Correlation Coefficients
#'
#' @description
#' Computes partial correlation coefficients between multiple features and a target
#' variable while controlling for an interference (confounding) variable. Adjusts
#' p-values for multiple testing using the Benjamini-Hochberg method.
#'
#' @param input Data frame containing feature variables, target variable, and
#'   interference variable.
#' @param interferenceid Character string specifying the column name of the
#'   interference (confounding) variable to control for.
#' @param target Character string specifying the column name of the target variable.
#' @param features Character vector specifying the column names of feature variables
#'   to correlate with the target.
#' @param method Character string specifying the correlation method. Options are
#'   \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}. Default is
#'   \code{"pearson"}.
#'
#' @return Tibble containing the following columns for each feature:
#' - `sig_names`: Feature name
#' - `p.value`: Raw p-value
#' - `statistic`: Partial correlation coefficient
#' - `p.adj`: Adjusted p-value (Benjamini-Hochberg method)
#' - `log10pvalue`: Negative log10-transformed p-value
#' - `stars`: Significance stars based on adjusted p-value thresholds
#'
#' @author Rongfang Shen
#' @export
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' sig_stad <- load_data("sig_stad")
#' # Calculate partial correlations controlling for tumor purity
#' res <- batch_pcc(
#'   input = sig_stad, interferenceid = "TumorPurity_estimate",
#'   target = "Pan_F_TBRs", method = "pearson",
#'   features = colnames(sig_stad)[70:ncol(sig_stad)]
#' )
batch_pcc <- function(input, interferenceid, target, features, method = "pearson") {
  # Input validation
  if (is.null(input) || !is.data.frame(input)) {
    stop("'input' must be a non-null data frame.")
  }
  if (!interferenceid %in% colnames(input)) {
    stop(sprintf("Interference variable '%s' not found in data columns.", interferenceid))
  }
  if (!target %in% colnames(input)) {
    stop(sprintf("Target variable '%s' not found in data columns.", target))
  }
  if (!is.character(features) || length(features) == 0) {
    stop("'features' must be a non-empty character vector.")
  }
  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop("'method' must be one of: 'pearson', 'spearman', 'kendall'.")
  }

  # Filter valid features
  features <- setdiff(features, c(interferenceid, target))
  features <- features[features %in% colnames(input)]
  if (length(features) == 0) {
    stop("None of the specified features were found in data columns.")
  }

  # Partial correlation test function
  pcor_test <- function(x, y, z, method = c("pearson", "spearman", "kendall")) {
    method <- match.arg(method)
    dat <- na.omit(cbind(x, y, z))
    if (nrow(dat) < 4) {
      return(c(estimate = NA_real_, p.value = NA_real_))
    }
    R <- cor(dat, method = method)
    rxy <- R["x", "y"]
    rxz <- R["x", "z"]
    ryz <- R["y", "z"]
    rho <- (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
    n <- nrow(dat)
    tstat <- rho * sqrt((n - 3) / (1 - rho^2))
    pval <- 2 * pt(abs(tstat), df = n - 3, lower.tail = FALSE)
    c(estimate = rho, p.value = pval)
  }

  # Vectorized partial correlation calculation
  results <- vapply(features, function(feat) {
    pcor_test(input[[feat]], input[[target]], input[[interferenceid]], method = method)
  }, numeric(2))

  # Build results data frame
  cc <- data.frame(
    sig_names = features,
    p.value = results["p.value", ],
    statistic = results["estimate", ],
    stringsAsFactors = FALSE
  )
  cc$p.adj <- p.adjust(cc$p.value, method = "BH")
  cc$log10pvalue <- -log10(cc$p.value)
  cc$stars <- cut(cc$p.adj,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
    labels = c("****", "***", "**", "*", "+", "")
  )
  cc <- cc[order(cc$p.value), , drop = FALSE]
  rownames(cc) <- NULL

  tibble::as_tibble(cc)
}
