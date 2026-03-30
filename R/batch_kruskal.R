#' Batch Kruskal-Wallis Test
#'
#' @description
#' Performs Kruskal-Wallis rank sum tests on multiple continuous features across
#' different groups. Computes p-values, adjusts for multiple testing, and ranks
#' features by significance.
#'
#' @param data Data frame containing the dataset for analysis.
#' @param group Character string specifying the name of the grouping variable.
#' @param feature Character vector specifying the names of feature variables to test.
#'   If \code{NULL}, the user is prompted to select all continuous features or specify
#'   features manually. Default is \code{NULL}.
#' @param feature_manipulation Logical indicating whether to apply feature manipulation
#'   to filter valid features. Default is \code{FALSE}.
#'
#' @return Tibble containing the following columns for each feature:
#' - `sig_names`: Feature name
#' - `p.value`: Raw p-value from Kruskal-Wallis test
#' - `p.adj`: Adjusted p-value (Benjamini-Hochberg method)
#' - `log10pvalue`: Negative log10-transformed p-value
#' - `stars`: Significance stars based on p-value thresholds
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' sig_stad <- load_data("sig_stad")
#' # Test features associated with TCGA molecular subtype
#' batch_kruskal(
#'   data = sig_stad, group = "Subtype",
#'   feature = colnames(sig_stad)[69:ncol(sig_stad)]
#' )
batch_kruskal <- function(data, group, feature = NULL, feature_manipulation = FALSE) {
  # Input validation
  if (is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a non-null data frame.")
  }
  if (nrow(data) == 0) {
    stop("'data' has no rows.")
  }
  if (!group %in% colnames(data)) {
    stop(sprintf("Group variable '%s' not found in data columns.", group))
  }

  # Prepare data
  data <- as.data.frame(data)
  data$group <- as.character(data[[group]])
  data <- data[!is.na(data$group) & data$group != "", , drop = FALSE]
  group_names <- unique(data$group)

  if (length(group_names) < 2) {
    stop("Group variable must have at least 2 unique values.")
  }

  # Feature selection
  if (is.null(feature)) {
    if (!interactive()) {
      stop("'feature' must be specified in non-interactive mode.")
    }
    message("'feature' must be specified, or all continuous features will be estimated...")
    index <- menu(c("all continuous features", "selected features"), title = "Choose features:")
    if (index == 1) {
      feature <- names(data)[vapply(data, is.numeric, logical(1))]
    } else {
      stop("Please specify the features that you want to proceed...")
    }
  }

  if (!is.character(feature) || length(feature) == 0) {
    stop("'feature' must be a non-empty character vector.")
  }

  # Filter valid features
  feature <- feature[feature %in% colnames(data)]
  if (length(feature) == 0) {
    stop("None of the specified features were found in data columns.")
  }

  if (feature_manipulation) {
    feature <- feature_manipulation(data = data, feature = feature, print_result = FALSE)
  }

  data <- data[, c("group", feature), drop = FALSE]

  message("Grouping information:")
  print(table(data$group))

  # Kruskal-Wallis test - vectorized
  results <- vapply(feature, function(feat) {
    test <- kruskal.test(data[[feat]] ~ data$group)
    c(statistic = test$statistic, p.value = test$p.value)
  }, numeric(2))

  cc <- data.frame(
    sig_names = feature,
    p.value = results["p.value", ],
    statistic = results["statistic", ],
    stringsAsFactors = FALSE
  )

  # Group means
  result_mean <- data %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise_if(is.numeric, mean, .groups = "drop")

  rownames(result_mean) <- result_mean$group
  result_mean$group <- NULL
  result_mean <- as.data.frame(t(result_mean))
  result_mean$sig_names <- rownames(result_mean)
  rownames(result_mean) <- NULL

  colnames(result_mean)[1:length(group_names)] <- sort(group_names)
  result_mean$mean <- rowMeans(result_mean[, sort(group_names), drop = FALSE], na.rm = TRUE)
  for (gn in sort(group_names)) {
    result_mean[[gn]] <- result_mean[[gn]] - result_mean$mean
  }

  # Merge and finalize
  cc <- merge(cc, result_mean, by = "sig_names", all.x = TRUE)
  cc$p.adj <- p.adjust(cc$p.value, method = "BH")
  cc$log10pvalue <- -log10(cc$p.value)
  cc$stars <- cut(cc$p.value,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
    labels = c("****", "***", "**", "*", "+", "")
  )
  cc <- cc[order(cc$p.value), , drop = FALSE]
  rownames(cc) <- NULL

  tibble::as_tibble(cc)
}
