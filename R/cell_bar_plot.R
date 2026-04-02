#' Visualize Cell Fractions as Stacked Bar Chart
#'
#' @description
#' Creates stacked bar charts to visualize tumor microenvironment (TME) cell
#' fractions. Supports batch visualization of deconvolution results from
#' methods such as CIBERSORT, EPIC, and quanTIseq.
#'
#' @param input Data frame containing deconvolution results.
#' @param id Character string specifying the column name containing sample
#'   identifiers. Default is "ID".
#' @param title Character string specifying the plot title.
#'   Default is "Cell Fraction".
#' @param features Character vector specifying column names representing cell
#'   types to plot. If NULL, columns are selected based on `pattern`.
#'   Default is NULL.
#' @param pattern Character string or regular expression to match column names
#'   for automatic feature selection. Used when `features` is NULL.
#'   Default is NULL.
#' @param legend.position Character string specifying legend position
#'   ("bottom", "top", "left", "right"). Default is "bottom".
#' @param coord_flip Logical indicating whether to flip plot coordinates using
#'   `coord_flip()`. Default is TRUE.
#' @param palette Integer specifying the color palette to use. Default is 3.
#' @param show_col Logical indicating whether to display color information.
#'   Default is FALSE.
#' @param cols Character vector of custom colors. If NULL, palette is used.
#'   Default is NULL.
#'
#' @return A ggplot2 object representing the stacked bar chart.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' set.seed(123)
#' input_data <- data.frame(
#'   ID = paste0("Sample", 1:10),
#'   Cell_A = runif(10, 0, 0.4),
#'   Cell_B = runif(10, 0, 0.3),
#'   Cell_C = runif(10, 0, 0.3)
#' )
#' cell_bar_plot(input = input_data, id = "ID", features = c("Cell_A", "Cell_B", "Cell_C"))
cell_bar_plot <- function(input, id = "ID", title = "Cell Fraction",
                          features = NULL, pattern = NULL,
                          legend.position = "bottom", coord_flip = TRUE,
                          palette = 3, show_col = FALSE, cols = NULL) {
  if (!is.data.frame(input)) {
    cli::cli_abort("{.arg input} must be a data frame")
  }
  if (nrow(input) == 0) {
    cli::cli_abort("{.arg input} has no rows")
  }
  if (!id %in% colnames(input)) {
    cli::cli_abort("ID column {.val {id}} not found in input")
  }

  input <- as.data.frame(input)
  colnames(input)[colnames(input) == id] <- "ID"

  if (is.null(features)) {
    if (is.null(pattern)) {
      cli::cli_abort("{.arg pattern} must be provided when {.arg features} is NULL")
    }
    feas <- colnames(input)[stringr::str_detect(colnames(input), pattern)]
    if (length(feas) == 0) {
      cli::cli_abort("No columns match pattern {.val {pattern}}")
    }
  } else {
    feas <- features[features %in% colnames(input)]
    if (length(feas) == 0) {
      cli::cli_abort("None of the specified features found in input")
    }
  }

  input <- input[, c("ID", feas), drop = FALSE]

  legend.direction <- if (legend.position %in% c("top", "bottom")) {
    "horizontal"
  } else {
    "vertical"
  }

  if (is.null(cols)) {
    cols <- palettes(
      category = "random",
      palette = palette %||% 4,
      show_col = show_col,
      show_message = TRUE
    )
  }

  plot_data <- tidyr::gather(input, "cell_type", "fraction", -.data$ID)

  lms <- rev(levels(input$ID))
  if (!all(lms %in% plot_data$ID)) {
    cli::cli_alert_warning("ID column is not a factor or has levels not present in data. Using unique IDs from data.")
    lms <- sort(unique(plot_data$ID))
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$ID, y = .data$fraction, fill = .data$cell_type
  )) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_light() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_x_discrete(limits = lms) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ggplot2::rel(2), hjust = 0.5),
      axis.text.x = ggplot2::element_text(
        face = "plain", angle = 0, hjust = 1, color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        face = "plain", angle = 30, hjust = 1, color = "black"
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.justification = c(0.5, 0.5),
      legend.box = "horizontal",
      legend.box.just = "top"
    )

  if (coord_flip) {
    p <- p + ggplot2::coord_flip()
  }

  p
}
