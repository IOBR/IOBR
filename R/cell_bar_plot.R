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
#' @param coord_filp Logical indicating whether to flip plot coordinates using
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
#' \donttest{
#' sig_stad <- load_data("sig_stad")
#' cell_bar_plot(
#'   input = sig_stad[1:20, ],
#'   id = "ID",
#'   features = colnames(sig_stad)[25:46]
#' )
#' }
cell_bar_plot <- function(input, id = "ID", title = "Cell Fraction",
                          features = NULL, pattern = NULL,
                          legend.position = "bottom", coord_flip = TRUE,
                          palette = 3, show_col = FALSE, cols = NULL) {
  input <- as.data.frame(input)
  colnames(input)[which(colnames(input) == id)] <- "ID"

  # Get features to plot
  if (is.null(features)) {
    if (is.null(pattern)) {
      cli::cli_abort("{.arg pattern} must be provided when {.arg features} is NULL")
    }
    feas <- colnames(input)[stringr::str_detect(colnames(input), pattern)]
  } else {
    feas <- features
  }

  input <- input[, c("ID", feas)]

  # Determine legend direction
  legend.direction <- if (legend.position %in% c("top", "bottom")) {
    "horizontal"
  } else {
    "vertical"
  }

  # Get colors
  if (is.null(cols)) {
    cols <- palettes(
      category = "random",
      palette = palette %||% 4,
      show_col = show_col,
      show_message = TRUE
    )
  }

  # Create plot data
  plot_data <- tidyr::gather(input, "cell_type", "fraction", -.data$ID)

  # Build plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$ID, y = .data$fraction, fill = .data$cell_type
  )) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_light() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_x_discrete(limits = rev(levels(input$ID))) +
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

  print(p)
  p
}
