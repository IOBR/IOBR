#' Calculate and Visualize Correlation Matrix Between Two Variable Sets
#'
#' @description
#' Calculates and visualizes the correlation matrix between two sets of variables.
#' Supports Pearson, Spearman, and Kendall correlation methods. The function
#' generates a customizable heatmap with significance stars.
#'
#' @param data Input data frame or matrix. Variables should be in columns.
#' @param feas1 Character vector of variable names for the first set.
#' @param feas2 Character vector of variable names for the second set.
#' @param method Correlation method: `"pearson"`, `"spearman"`, or `"kendall"`.
#'   Default is `"pearson"`.
#' @param path Directory to save the plot. If `NULL`, plot is not saved.
#'   Default is `NULL`.
#' @param index Numeric prefix for output filename. Default is 1.
#' @param fig.type File format: `"pdf"`, `"png"`, etc. Default is `"pdf"`.
#' @param width Plot width in inches. Auto-calculated if `NULL`.
#' @param height Plot height in inches. Auto-calculated if `NULL`.
#' @param project Project name for plot title. Default is `NULL`.
#' @param is.matrix Logical: if `TRUE`, data is transposed. Default is `FALSE`.
#' @param scale Logical: scale variables before correlation. Default is `TRUE`.
#' @param font.size Font size for axis labels. Default is 15.
#' @param fill_by_cor Logical: show correlation values instead of stars.
#'   Default is `FALSE`.
#' @param round.num Decimal places for correlation values. Default is 1.
#' @param font.size.star Font size for significance stars. Default is 8.
#' @param cols Custom colors for gradient (low, mid, high). If `NULL`, uses
#'   blue-white-red. Default is `NULL`.
#'
#' @return ggplot object displaying the correlation matrix heatmap.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' data <- as.data.frame(matrix(rnorm(1000), nrow = 100, ncol = 10))
#' colnames(data) <- paste0("Gene_", 1:10)
#'
#' feas1 <- c("Gene_1", "Gene_2", "Gene_3")
#' feas2 <- c("Gene_4", "Gene_5", "Gene_6")
#'
#' cor_plot <- get_cor_matrix(
#'   data = data,
#'   feas1 = feas1,
#'   feas2 = feas2,
#'   method = "spearman",
#'   project = "Example Correlation"
#' )
#' }
get_cor_matrix <- function(data,
                           feas1,
                           feas2,
                           method = c("pearson", "spearman", "kendall"),
                           path = NULL,
                           index = 1,
                           fig.type = "pdf",
                           width = NULL,
                           height = NULL,
                           project = NULL,
                           is.matrix = FALSE,
                           scale = TRUE,
                           font.size = 15,
                           fill_by_cor = FALSE,
                           round.num = 1,
                           font.size.star = 8,
                           cols = NULL) {
  # Validate method
  method <- rlang::arg_match(method)

  # Create output folder if requested
  if (!is.null(path)) {
    path <- creat_folder(path)
  }

  # Handle matrix input
  if (is.matrix) data <- as.data.frame(t(data))

  # Validate features
  if (!is.character(feas1) || length(feas1) == 0) {
    cli::cli_abort("{.arg feas1} must be a non-empty character vector.")
  }
  if (!is.character(feas2) || length(feas2) == 0) {
    cli::cli_abort("{.arg feas2} must be a non-empty character vector.")
  }

  # Filter to existing features
  feas1 <- feas1[feas1 %in% colnames(data)]
  feas2 <- feas2[feas2 %in% colnames(data)]

  if (length(feas1) == 0) {
    cli::cli_abort("No features from {.arg feas1} found in data.")
  }
  if (length(feas2) == 0) {
    cli::cli_abort("No features from {.arg feas2} found in data.")
  }

  cli::cli_alert_info(
    "Calculating {method} correlation: {length(feas1)} x {length(feas2)}"
  )

  # Scale data if requested
  data_vars <- data[, unique(c(feas1, feas2)), drop = FALSE]
  if (scale) {
    data_vars <- scale(data_vars)
  }

  # Calculate correlation
  rlang::check_installed("psych", reason = "to calculate correlation matrices")
  result <- psych::corr.test(
    data_vars[, feas1, drop = FALSE],
    data_vars[, feas2, drop = FALSE],
    method = method
  )

  # Reshape for plotting
  heat <- cbind(
    reshape2::melt(result$r),
    reshape2::melt(result$p)
  )[, c(1, 2, 3, 6)]
  colnames(heat) <- c("ID1", "ID2", "cor", "pvalue")

  # Clean labels for display (replace underscores with spaces)
  feas1_clean <- gsub("_", " ", feas1)
  feas2_clean <- gsub("_", " ", feas2)

  heat$ID1 <- gsub("_", " ", as.character(heat$ID1))
  heat$ID2 <- gsub("_", " ", as.character(heat$ID2))

  # Set factor levels to preserve order
  heat$ID1 <- factor(heat$ID1, levels = feas1_clean)
  heat$ID2 <- factor(heat$ID2, levels = rev(feas2_clean))

  # Create labels
  if (fill_by_cor) {
    heat$stars <- round(heat$cor, round.num)
  } else {
    heat$stars <- cut(
      heat$pvalue,
      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.5, Inf),
      labels = c("***", "**", "*", "+", "")
    )
  }

  # Define colors
  if (is.null(cols)) {
    low_col <- "#2C7BB6" # Blue
    mid_col <- "white"
    high_col <- "#D7191C" # Red
  } else {
    if (length(cols) < 2) {
      cli::cli_abort("{.arg cols} must have at least 2 colors (low, high)")
    }
    low_col <- cols[1]
    high_col <- cols[length(cols)]
    mid_col <- if (length(cols) >= 3) cols[2] else "white"
  }

  # Create plot
  p <- ggplot2::ggplot(heat, ggplot2::aes(x = .data$ID1, y = .data$ID2, fill = .data$cor))

  cor_plot <- p +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = low_col,
      mid = mid_col,
      high = high_col,
      name = "Coefficient"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$stars),
      color = "black",
      size = font.size.star
    ) +
    ggplot2::scale_y_discrete(
      labels = function(y) stringr::str_wrap(y, width = 40)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = -45,
        hjust = 0,
        size = font.size
      ),
      axis.text.y = ggplot2::element_text(
        angle = 0,
        hjust = 1,
        size = font.size
      ),
      axis.title = ggplot2::element_text(
        size = font.size + 4
      )
    )

  if (!is.null(project)) {
    cor_plot <- cor_plot + ggplot2::ggtitle(label = project)
  }

  # Auto-calculate dimensions
  if (is.null(width)) width <- length(feas1) * 0.55 + 6.5
  if (is.null(height)) height <- length(feas2) * 0.35 + 4.5

  # Save if path provided
  if (!is.null(path)) {
    ggplot2::ggsave(
      cor_plot,
      filename = paste0(index, "-", project %||% "cor", "-cor_plot.", fig.type),
      width = width,
      height = height,
      path = path$folder_name
    )
  }

  print(cor_plot)
  invisible(cor_plot)
}
