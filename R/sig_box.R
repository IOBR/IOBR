#' Signature Box Plot with Statistical Comparisons
#'
#' @description
#' Creates box plots to visualize signature distributions across groups with
#' optional statistical pairwise comparisons. Supports both data frames and
#' Seurat objects for single-cell data visualization.
#'
#' @param data Data frame or Seurat object containing the signature and
#'   grouping variable.
#' @param signature Character string specifying the column name (or feature name
#'   in Seurat) for the signature values to plot on the y-axis.
#' @param variable Character string specifying the grouping variable column name
#'   for the x-axis.
#' @param palette Character string specifying the color palette name.
#'   Default is `"nrc"`.
#' @param cols Character vector of custom fill colors. If `NULL`, palette is
#'   used. Default is `NULL`.
#' @param jitter Logical indicating whether to add jittered points to the box
#'   plot. Default is `FALSE`.
#' @param point_size Numeric value specifying the size of jittered points.
#'   Default is 5.
#' @param angle_x_text Numeric value specifying the rotation angle for x-axis
#'   labels (in degrees). Default is 0.
#' @param hjust Numeric value specifying the horizontal justification of x-axis
#'   labels. Default is 0.5.
#' @param show_pairwise_p Logical indicating whether to display pairwise
#'   comparison p-values between groups. Default is `TRUE`.
#' @param show_overall_p Logical indicating whether to display the overall
#'   group difference p-value. Default is `FALSE`.
#' @param return_stat_res Logical indicating whether to return statistical test
#'   results instead of the plot. Default is `FALSE`.
#' @param size_of_pvalue Numeric value specifying the font size for p-values.
#'   Default is 6.
#' @param size_of_font Numeric value specifying the base font size.
#'   Default is 10.
#' @param assay Character string specifying the assay name (for Seurat objects).
#'   Default is `NULL`.
#' @param slot Character string specifying the slot name (for Seurat objects).
#'   Default is `"scale.data"`.
#' @param scale Logical indicating whether to scale signature values (z-score
#'   transformation). Default is `FALSE`.
#'
#' @return If `return_stat_res = FALSE`, returns a ggplot2 object. If
#'   `return_stat_res = TRUE`, returns a data frame containing statistical test
#'   results.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' # Create small example data
#' set.seed(123)
#' test_data <- data.frame(
#'   subtype = rep(c("A", "B"), each = 50),
#'   TMEscore_plus = rnorm(100)
#' )
#' sig_box(
#'   data = test_data,
#'   signature = "TMEscore_plus",
#'   variable = "subtype",
#'   jitter = TRUE,
#'   palette = "jco"
#' )
sig_box <- function(data,
                    signature,
                    variable,
                    palette = "nrc",
                    cols = NULL,
                    jitter = FALSE,
                    point_size = 5,
                    angle_x_text = 0,
                    hjust = 0.5,
                    show_pairwise_p = TRUE,
                    show_overall_p = FALSE,
                    return_stat_res = FALSE,
                    size_of_pvalue = 6,
                    size_of_font = 10,
                    assay = NULL,
                    slot = "scale.data",
                    scale = FALSE) {
  if (is.null(data)) return(NULL)
  rlang::check_installed("ggpubr")

  # Handle Seurat object
  if (inherits(data, "Seurat")) {
    cli::cli_alert_info("Extracting data from Seurat object...")
    data <- extract_sc_data(
      sce = data,
      vars = signature,
      assay = assay,
      slot = slot,
      combine_meta_data = TRUE
    )
  }

  # Input validation
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame or Seurat object")
  }
  if (!signature %in% colnames(data)) {
    cli::cli_abort("Signature {.val {signature}} not found in data")
  }
  if (!variable %in% colnames(data)) {
    cli::cli_abort("Variable {.val {variable}} not found in data")
  }

  data <- as.data.frame(data)
  data <- data[, c(variable, signature), drop = FALSE]
  colnames(data)[colnames(data) == variable] <- "variable"
  colnames(data)[colnames(data) == signature] <- "signature"

  # Scale if requested
  if (scale) {
    data$signature <- as.numeric(scale(data$signature))
  }

  # Remove NA in grouping variable or signature
  data <- data[!is.na(data$variable) & !is.na(data$signature), , drop = FALSE]

  if (nrow(data) == 0) {
    cli::cli_abort("No valid data after removing NA values")
  }

  # Get colors
  if (is.null(cols)) {
    cols <- palettes(
      category = "box",
      palette = palette,
      show_message = FALSE,
      show_col = FALSE
    )
  }

  # Create base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data$variable,
    y = .data$signature,
    fill = .data$variable
  )) +
    ggplot2::geom_boxplot(
      notch = FALSE,
      outlier.shape = NA
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::ylab(signature) +
    ggplot2::xlab(variable)

    # Count samples per group
  group_count <- table(data$variable)
  valid_groups <- names(group_count[group_count >= 2])

  # Generate comparisons safely
  if (length(valid_groups) >= 2) {
    comparisons <- utils::combn(valid_groups, 2, simplify = FALSE)
  } else {
    comparisons <- list()
  }

  # Adjust hjust for angled text
  if (angle_x_text %in% c(30, 45, 60)) {
    hjust <- 1
  }

  # Apply theme
  size_font <- size_of_font * 0.2
  p <- p +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(size_font)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(size_font)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(size_font)),
      axis.text.x = ggplot2::element_text(
        face = "plain",
        angle = angle_x_text,
        hjust = hjust,
        color = "black"
      ),
      axis.line = ggplot2::element_line(color = "grey", size = 0.05),
      legend.position = "none"
    )

  # Add statistical comparisons
  if (show_pairwise_p) {
    if (length(comparisons) > 0) {
      p <- p + ggpubr::stat_compare_means(
        comparisons = comparisons,
        size = size_of_pvalue
      )
    }
  }

  if (show_overall_p) {
    if (length(valid_groups) >= 2) {
      stat_data <- data[data$variable %in% valid_groups, , drop = FALSE]
      p <- p + ggpubr::stat_compare_means(
        data = stat_data,
        size = size_of_pvalue
      )
    }
  }

  # Calculate and print statistics
  res <- NULL
  if (length(valid_groups) >= 2) {
    stat_data <- data[data$variable %in% valid_groups, , drop = FALSE]

    res <- tryCatch({
      ggpubr::compare_means(
        stats::reformulate("variable", "signature"),
        data = stat_data
      )
    }, error = function(e) {
      NULL
    })

    if (interactive() && !is.null(res)) print(res)
  }
  
   # Add jitter points if requested
  if (jitter) {
    p <- p + ggplot2::geom_point(
      shape = 21,
      size = point_size,
      position = ggplot2::position_jitterdodge(dodge.width = 0.2),
      alpha = 0.5
    )
  }

  if (return_stat_res) {
    invisible(res)
  } else {
    invisible(p)
  }
}
