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
#' \itemize{
#'   \item \code{sig_names}: Feature name
#'   \item \code{p.value}: Raw p-value
#'   \item \code{statistic}: Partial correlation coefficient
#'   \item \code{p.adj}: Adjusted p-value (Benjamini-Hochberg method)
#'   \item \code{log10pvalue}: Negative log10-transformed p-value
#'   \item \code{stars}: Significance stars based on adjusted p-value thresholds
#' }
#'
#' @author Rongfang Shen
#' @export
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Calculate partial correlations controlling for tumor purity
#' res <- batch_pcc(input = sig_stad, interferenceid = "TumorPurity_estimate",
#'                  target = "Pan_F_TBRs", method = "pearson",
#'                  features = colnames(sig_stad)[70:ncol(sig_stad)])
batch_pcc <- function(input, interferenceid, target, features, method = "pearson") {
  dat <- input
  features <- setdiff(features, c(unique(interferenceid, target)))

  aa <- dat[, features] %>%
    tibble::as_tibble() %>%
    map(ppcor::pcor.test, y = dat[, target], z = dat[, interferenceid], method = method)
  pvalue <- aa %>% purrr::map_dbl("p.value")
  statistic <- aa %>% purrr::map_dbl("estimate")
  cc <- data.frame(
    sig_names = features,
    p.value = pvalue,
    statistic = statistic
  )
  cc <- cc[order(cc$p.value, decreasing = F), ]
  cc$p.adj <- p.adjust(cc$p.value, method = "BH")
  cc$log10pvalue <- -1 * log10(cc$p.value)
  rownames(cc) <- NULL
  cc$stars <- cut(cc$p.adj,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
    label = c("****", "***", "**", "*", "+", "")
  )
  cc <- tibble::as_tibble(cc)
  return(cc)
}
