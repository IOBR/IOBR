#' Batch Correlation Analysis
#'
#' @description
#' Performs correlation analysis between a target variable and multiple feature
#' variables. Computes correlation coefficients, p-values, and adjusts for multiple
#' testing using the Benjamini-Hochberg method.
#'
#' @param data Data frame containing the target and feature variables.
#' @param target Character string specifying the name of the target variable.
#' @param feature Character vector specifying the names of feature variables to
#'   correlate with the target.
#' @param method Character string specifying the correlation method. Options are
#'   \code{"spearman"}, \code{"pearson"}, or \code{"kendall"}. Default is
#'   \code{"spearman"}.
#'
#' @return Tibble containing the following columns for each feature:
#' \itemize{
#'   \item \code{sig_names}: Feature name
#'   \item \code{p.value}: Raw p-value
#'   \item \code{statistic}: Correlation coefficient
#'   \item \code{p.adj}: Adjusted p-value (Benjamini-Hochberg method)
#'   \item \code{log10pvalue}: Negative log10-transformed p-value
#'   \item \code{stars}: Significance stars based on p-value thresholds
#' }
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Perform batch correlation analysis
#' results <- batch_cor(data = sig_stad, target = "CD_8_T_effector",
#'                      feature = colnames(sig_stad)[69:ncol(sig_stad)])
batch_cor <- function(data, target, feature, method = "spearman") {
  if (!target %in% colnames(data)) stop(">>>== target was not in the colnames of data")
  data <- as.data.frame(data)
  feature <- feature[feature %in% colnames(data)]

  if (length(feature) == 0) stop(">>>== features were not in the colnames of data")

  feature <- feature_manipulation(data = data, feature = feature)

  feature <- feature[!feature == target]

  data <- data[!is.na(data[, target]), ]

  feature <- feature[sapply(data[, feature], function(x) sd(x, na.rm = TRUE) > 0)]

  aa <- lapply(data[, feature], function(x) cor.test(x, data[, target], method = method))

  bb <- lapply(data[, feature], function(x) exact_pvalue(x, data[, target], method = method))

  bb <- as.data.frame(bb)

  bb <- as.data.frame(t(bb))
  # print(head(bb))
  cc <- data.frame(
    sig_names = feature,
    p.value = bb$V1,
    statistic = sapply(aa, getElement, name = "estimate")
  )

  # cc<-cc[base::order(cc$p.value, decreasing = FALSE),]
  cc$p.adj <- p.adjust(cc$p.value, method = "BH")
  cc$log10pvalue <- -1 * log10(cc$p.value)
  rownames(cc) <- NULL
  cc$stars <- cut(cc$p.value,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
    label = c("****", "***", "**", "*", "+", "")
  )
  cc <- cc[base::order(cc$p.value, decreasing = FALSE), ]
  cc <- tibble::as_tibble(cc)
  print(head(cc))
  return(cc)
}


#' Calculate Exact P-Value for Correlation
#'
#' @description
#' Computes the exact p-value for the correlation between two numeric variables
#' using a specified correlation method. This function provides detailed statistical
#' support for correlation analyses.
#'
#' @param x Numeric vector representing the first variable.
#' @param y Numeric vector representing the second variable.
#' @param method Character string specifying the correlation method. Options include
#'   \code{"spearman"}, \code{"pearson"}, or \code{"kendall"}.
#'
#' @return Numeric value representing the exact p-value from the correlation test.
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Calculate exact p-value for correlation between two variables
#' p_val <- exact_pvalue(x = sig_stad$CD8.T.cells, y = sig_stad$CD_8_T_effector,
#'                       method = "spearman")
#' print(p_val)
exact_pvalue <- function(x, y, method) {
  l <- length(y)
  r <- cor(x = x, y = y, method = method, use = "complete.obs")

  if (r < 0) {
    t <- r * sqrt((l - 2) / (1 - r^2))
    s <- (l^3 - l) * (1 - r) / 6
    p <- 2 * (pt(q = t, df = l - 2))
  } else {
    t <- r * sqrt((l - 2) / (1 - r^2))
    s <- (l^3 - l) * (1 - r) / 6
    p <- 2 * (pt(q = -t, df = l - 2))
  }
  return(p)
}
