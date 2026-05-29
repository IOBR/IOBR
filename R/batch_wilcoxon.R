#' Batch Wilcoxon Rank-Sum Test Between Two Groups
#'
#' @description
#' Performs Wilcoxon rank-sum tests (Mann-Whitney U tests) to compare the
#' distribution of specified features between two groups. Computes p-values,
#' adjusts for multiple testing, and ranks features by significance.
#'
#' @param data Data frame containing the dataset for analysis.
#' @param target Character string specifying the column name of the grouping
#'   variable. Default is `"group"`.
#' @param feature Character vector specifying feature names to analyze. If
#'   `NULL`, prompts for selection (interactive mode only). Default is `NULL`.
#' @param feature_manipulation Logical indicating whether to apply feature
#'   manipulation filtering. Default is `FALSE`.
#'
#' @return Tibble with columns:
#' \describe{
#'   \item{sig_names}{Feature name}
#'   \item{p.value}{Raw p-value}
#'   \item{statistic}{Difference in means between groups}
#'   \item{p.adj}{Adjusted p-value (Benjamini-Hochberg)}
#'   \item{log10pvalue}{Negative log10-transformed p-value}
#'   \item{stars}{Significance stars: **** p<0.0001, *** p<0.001,
#'     ** p<0.01, * p<0.05, + p<0.5}
#'   \item{group1, group2}{Mean values for each group}
#' }
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' # Create small example data
#' set.seed(123)
#' test_data <- data.frame(
#'   Gender = rep(c("Male", "Female"), each = 50),
#'   Signature1 = rnorm(100),
#'   Signature2 = rnorm(100)
#' )
#' # Compare features by gender
#' res <- batch_wilcoxon(
#'   data = test_data,
#'   target = "Gender",
#'   feature = c("Signature1", "Signature2")
#' )
#' head(res)
batch_wilcoxon <- function(data,
                           target = "group",
                           feature = NULL,
                           feature_manipulation = FALSE) {
  if (is.null(data)) return(NULL)
  # Input validation
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a non-null data frame.")
  }
  if (!target %in% colnames(data)) {
    cli::cli_abort("Target variable {.val {target}} not found in data.")
  }

  data <- as.data.frame(data)

  # Rename target column
  colnames(data)[colnames(data) == target] <- "group"

  # Clean group data
  data <- data[!is.na(data$group) & data$group != "", , drop = FALSE]
  data$group <- as.character(data$group)

  group_names <- unique(data$group)
  if (length(group_names) != 2) {
    cli::cli_abort(c(
      "Wilcoxon test requires exactly 2 groups.",
      "i" = "Found {length(group_names)} group{?s}: {.val {group_names}}"
    ))
  }

  # Feature selection
  if (is.null(feature)) {
    if (!interactive()) {
      cli::cli_abort(paste(
        "{.arg feature} must be specified in",
        "non-interactive mode."
      ))
    }
    cli::cli_alert_info("Select features for analysis")
    index <- utils::menu(
      c("All continuous features", "Specify manually"),
      title = "Feature selection:"
    )
    if (index == 1) {
      feature <- colnames(data)[vapply(data, is.numeric, logical(1))]
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
    cli::cli_alert_warning(paste(
      "Ignoring {length(invalid_features)} invalid feature{?s}:",
      "{.val {invalid_features}}"
    ))
  }

  if (length(valid_features) == 0) {
    cli::cli_abort("No valid features found in data.")
  }

  # Apply feature manipulation if requested
  if (feature_manipulation) {
    valid_features <- IOBR::feature_manipulation(
      data = data,
      feature = valid_features,
      print_result = FALSE
    )
  }

  cli::cli_alert_info("Groups: {.val {group_names}}")
  cli::cli_alert_info("Features: {length(valid_features)}")

  # Subset data
  data <- data[, c("group", valid_features), drop = FALSE]

  # Perform Wilcoxon tests
  test_results <- lapply(valid_features, function(feat) {
    stats::wilcox.test(data[[feat]] ~ data$group, exact = FALSE)
  })
  names(test_results) <- valid_features

  # Calculate group means - vectorized approach
  result_mean <- data %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(valid_features), function(.x) mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  # Convert to wide format with groups as columns
  result_mean <- tidyr::pivot_longer(result_mean,
    cols = valid_features
  ) %>% tidyr::pivot_wider(
    id_cols = "name", names_from = "group"
  )

  # Calculate statistic (difference between groups)
  if (ncol(result_mean) >= 3) { # name + 2 groups column
    result_mean$statistic <- result_mean[[2]] - result_mean[[3]]
  } else {
    result_mean$statistic <- NA_real_
  }

  # Build results
  cc <- data.frame(
    name = valid_features,
    p.value = vapply(test_results, function(x) x$p.value, numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  cc <- dplyr::left_join(cc, result_mean, by = "name") %>%
    dplyr::arrange(.data$p.value) %>%
    dplyr::mutate(
      p.adj = stats::p.adjust(.data$p.value, method = "BH"),
      log10pvalue = -log10(.data$p.value),
      stars = cut(.data$p.value,
        breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
        labels = c("****", "***", "**", "*", "+", "")
      )
    )
  colnames(cc)[1] <- "sig_names"

  cli::cli_alert_success("Wilcoxon test complete")

  tibble::as_tibble(cc)
}
