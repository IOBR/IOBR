#' Forest Plot for Survival Analysis Results
#'
#' @description
#' Generates a forest plot to visualize hazard ratios, confidence intervals,
#' and p-values for gene signatures or features from survival analysis.
#'
#' @param data Data frame with survival analysis results including p-values,
#'   hazard ratios, and confidence intervals.
#' @param signature Character string. Column name for signatures or feature names.
#' @param pvalue Character string. Column name for p-values. Default is `"P"`.
#' @param HR Character string. Column name for hazard ratios. Default is `"HR"`.
#' @param CI_low_0.95 Character string. Column name for lower CI bound.
#'   Default is `"CI_low_0.95"`.
#' @param CI_up_0.95 Character string. Column name for upper CI bound.
#'   Default is `"CI_up_0.95"`.
#' @param n Integer. Maximum number of signatures to display. Default is `10`.
#' @param max_character Integer. Maximum characters for labels before wrapping.
#'   Default is `25`.
#' @param discrete_width Integer. Width for discretizing long labels. Default is `35`.
#' @param color_option Integer. Color option for p-value gradient (1, 2, or 3).
#'   Default is `1`.
#' @param cols Character vector. Custom colors for p-value gradient (low to high).
#'   Default is `NULL`.
#' @param text.size Numeric. Text size for y-axis labels. Default is `13`.
#'
#' @return A ggplot2 object of the forest plot.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' # Example with sample survival results
#' sample_results <- data.frame(
#'   ID = c("Sig1", "Sig2", "Sig3"),
#'   HR = c(1.5, 0.8, 2.0),
#'   P = c(0.01, 0.05, 0.001),
#'   CI_low_0.95 = c(1.1, 0.6, 1.5),
#'   CI_up_0.95 = c(2.0, 1.0, 2.8)
#' )
#' sig_forest(data = sample_results, signature = "ID")
sig_forest <- function(data,
                       signature,
                       pvalue = "P",
                       HR = "HR",
                       CI_low_0.95 = "CI_low_0.95",
                       CI_up_0.95 = "CI_up_0.95",
                       n = 10,
                       max_character = 25,
                       discrete_width = 35,
                       color_option = 1,
                       cols = NULL,
                       text.size = 13) {
  rlang::check_installed("stringr")

  data <- as.data.frame(data)

  # Validate required columns
  required_cols <- c(signature, pvalue, HR, CI_low_0.95, CI_up_0.95)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing required columns: {.val {missing_cols}}")
  }

  # Standardize column names
  colnames(data)[colnames(data) == signature] <- "signature"
  colnames(data)[colnames(data) == pvalue] <- "P"
  colnames(data)[colnames(data) == HR] <- "HR"
  colnames(data)[colnames(data) == CI_low_0.95] <- "CI_low_0.95"
  colnames(data)[colnames(data) == CI_up_0.95] <- "CI_up_0.95"

  # Convert to numeric and remove incomplete cases
  data[, c("P", "HR", "CI_low_0.95", "CI_up_0.95")] <- lapply(
    data[, c("P", "HR", "CI_low_0.95", "CI_up_0.95")],
    as.numeric
  )
  data <- data[stats::complete.cases(data), ]

  if (nrow(data) == 0) {
    cli::cli_abort("No valid data after removing incomplete cases")
  }

  # Select top signatures
  if (nrow(data) > n) {
    cli::cli_alert_info("Showing top {n} signatures")

    data$HR_statistic <- data$HR - 1
    good_features <- high_var_fea(
      result = data,
      target = "signature",
      name_padj = "P",
      padj_cutoff = 1,
      name_logfc = "HR_statistic",
      logfc_cutoff = 0,
      n = n / 2
    )
    data <- data[data$signature %in% good_features, ]
    data <- data[order(data$HR, decreasing = FALSE), ]
    data$signature <- as.character(data$signature)
  }

  # Process long labels
  if (max(nchar(as.character(data$signature))) > max_character) {
    long_idx <- nchar(as.character(data$signature)) > max_character
    data$signature[long_idx] <- gsub("_", " ", as.character(data$signature[long_idx]))
  }

  # Set factor order based on HR
  data$signature <- factor(data$signature, levels = data$signature)

  # Set gradient colors
  if (!is.null(cols)) {
    if (length(cols) < 2) {
      cli::cli_abort("{.arg cols} must have at least 2 colors for gradient")
    }
    gradient_colors <- grDevices::colorRampPalette(cols)(256)
  } else {
    gradient_colors <- grDevices::colorRampPalette(c(
      "#000004FF", "#51127CFF", "#B63679FF", "#FCA50AFF", "#F7F419FF"
    ))(256)
  }

  # Create forest plot
  pp <- ggplot2::ggplot(
    data = data,
    ggplot2::aes(x = .data$HR, y = .data$signature, color = .data$P)
  ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmax = .data$CI_up_0.95, xmin = .data$CI_low_0.95),
      color = "black",
      height = 0,
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$HR, y = .data$signature),
      size = 4.5,
      shape = 16
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::scale_x_continuous(breaks = c(0.5, 1, 1.50)) +
    ggplot2::coord_trans(x = "log2") +
    ggplot2::ylab("Features") +
    ggplot2::xlab("Hazard Ratios of Features") +
    ggplot2::labs(color = "P value") +
    ggplot2::scale_color_gradientn(colors = gradient_colors, name = "P value") +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 15, color = "black", vjust = 0.5, hjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        size = text.size, color = "black", vjust = 0.5, hjust = 1
      ),
      plot.title = ggplot2::element_text(
        size = 15, colour = "black", vjust = 0.5, hjust = 0.5
      )
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x) stringr::str_wrap(x, width = discrete_width)
    )

  pp
}
