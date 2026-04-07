#' Batch Kruskal-Wallis Test
#'
#' @description
#' Performs Kruskal-Wallis rank sum tests on multiple continuous features across
#' different groups. Computes p-values, adjusts for multiple testing, and ranks
#' features by significance.
#'
#' @param data Data frame containing the dataset for analysis.
#' @param group Character string specifying the name of the grouping variable.
#' @param feature Character vector specifying the names of feature variables to
#'   test. If `NULL`, the user is prompted to select features (interactive mode
#'   only). Default is `NULL`.
#' @param feature_manipulation Logical indicating whether to apply feature
#'   manipulation to filter valid features. Default is `FALSE`.
#'
#' @return Tibble containing:
#' \describe{
#'   \item{sig_names}{Feature name}
#'   \item{p.value}{Raw p-value from Kruskal-Wallis test}
#'   \item{statistic}{Test statistic (chi-squared)}
#'   \item{p.adj}{Adjusted p-value (Benjamini-Hochberg)}
#'   \item{log10pvalue}{Negative log10-transformed p-value}
#'   \item{stars}{Significance stars: **** p<0.0001, *** p<0.001,
#'     ** p<0.01, * p<0.05, + p<0.5}
#'   \item{group columns}{Mean-centered values for each group}
#' }
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' # Load TCGA-STAD signature data
#' sig_stad <- load_data("sig_stad")
#'
#' # Test features by gender (if available in your dataset)
#' if ("Gender" %in% colnames(sig_stad)) {
#'   res <- batch_kruskal(
#'     data = sig_stad,
#'     group = "Gender",
#'     feature = colnames(sig_stad)[69:ncol(sig_stad)]
#'   )
#'   head(res)
#' }
batch_kruskal <- function(data,
                          group,
                          feature = NULL,
                          feature_manipulation = FALSE) {
  # Input validation
  if (is.null(data) || !is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a non-null data frame.")
  }
  if (nrow(data) == 0) {
    cli::cli_abort("{.arg data} has no rows.")
  }
  if (!group %in% colnames(data)) {
    cli::cli_abort("Group variable {.val {group}} not found in data.")
  }

  # Prepare data
  data <- as.data.frame(data)
  data$group <- as.character(data[[group]])
  data <- data[!is.na(data$group) & data$group != "", , drop = FALSE]

  group_names <- unique(data$group)
  if (length(group_names) < 2) {
    cli::cli_abort("Group variable must have at least 2 unique non-NA values.")
  }

  # Feature selection
  if (is.null(feature)) {
    if (!interactive()) {
      cli::cli_abort(
        "{.arg feature} must be specified in non-interactive mode."
      )
    }
    cli::cli_alert_info("Select features for analysis")
    index <- utils::menu(
      c("All continuous features", "Specify manually"),
      title = "Feature selection:"
    )
    if (index == 1) {
      feature <- names(data)[vapply(data, is.numeric, logical(1))]
    } else {
      cli::cli_abort("Please specify features via the {.arg feature} parameter")
    }
  }

  if (!is.character(feature) || length(feature) == 0) {
    cli::cli_abort("{.arg feature} must be a non-empty character vector.")
  }

  # Filter valid features
  valid_features <- feature[feature %in% colnames(data)]
  invalid_features <- setdiff(feature, colnames(data))

  if (length(invalid_features) > 0) {
    cli::cli_alert_warning(
      "Ignoring {length(invalid_features)} invalid feature{?s}:",
      " {.val {invalid_features}}"
    )
  }

  if (length(valid_features) == 0) {
    cli::cli_abort("No valid features found in data.")
  }

  # Apply feature manipulation if requested
  if (feature_manipulation) {
    valid_features <- feature_manipulation(
      data = data,
      feature = valid_features,
      print_result = FALSE
    )
  }

  data <- data[, c("group", valid_features), drop = FALSE]

  cli::cli_alert_info("Groups: {length(group_names)} ({.val {group_names}})")
  cli::cli_alert_info("Features: {length(valid_features)}")

  # Kruskal-Wallis test - vectorized
  results <- vapply(valid_features, function(feat) {
    test <- stats::kruskal.test(data[[feat]] ~ data$group)
    c(statistic = as.numeric(test$statistic), p.value = test$p.value)
  }, numeric(2))

  cc <- data.frame(
    sig_names = valid_features,
    p.value = results["p.value", ],
    statistic = results["statistic", ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Calculate group means (mean-centered)
  result_mean <- data %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(valid_features),
        function(.x) mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(valid_features),
      names_to = "sig_names",
      values_to = "value"
    ) %>%
    tidyr::pivot_wider(
      names_from = "group",
      values_from = "value"
    ) %>%
    as.data.frame()

  # Calculate mean-centered values
  group_cols <- intersect(group_names, colnames(result_mean))
  if (length(group_cols) > 0) {
    result_mean$overall_mean <- rowMeans(
      as.matrix(result_mean[, group_cols, drop = FALSE]),
      na.rm = TRUE
    )

    for (gn in group_cols) {
      result_mean[[gn]] <- result_mean[[gn]] - result_mean$overall_mean
    }
    result_mean$overall_mean <- NULL
  }

  # Merge results
  cc <- merge(cc, result_mean, by = "sig_names", all.x = TRUE)
  cc$p.adj <- stats::p.adjust(cc$p.value, method = "BH")
  cc$log10pvalue <- -log10(cc$p.value)
  cc$stars <- cut(cc$p.value,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
    labels = c("****", "***", "**", "*", "+", "")
  )

  cc <- cc[order(cc$p.value), , drop = FALSE]

  cli::cli_alert_success("Kruskal-Wallis test complete")

  tibble::as_tibble(cc)
}
