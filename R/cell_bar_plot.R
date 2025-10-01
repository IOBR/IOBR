#' Visualize Cell Fractions as Stacked Bar Chart
#'
#' @description
#' Creates stacked bar charts to visualize tumor microenvironment (TME) cell fractions.
#' Supports batch visualization of deconvolution results from methods such as CIBERSORT,
#' EPIC, and quanTIseq. Enables comparison of TME cell distributions across samples.
#'
#' @param input Data frame containing deconvolution results (e.g., from CIBERSORT,
#'   quanTIseq, or EPIC).
#' @param id Character string specifying the column name containing sample identifiers.
#'   Default is \code{"ID"}.
#' @param title Character string specifying the plot title. Default is \code{"Cell Fraction"}.
#' @param features Character vector specifying column names representing cell types to
#'   plot. If \code{NULL}, columns are selected based on \code{pattern}. Default is
#'   \code{NULL}.
#' @param pattern Character string or regular expression to match column names for
#'   automatic feature selection. Used when \code{features} is \code{NULL}. Default is
#'   \code{NULL}.
#' @param legend.position Character string specifying legend position (\code{"bottom"},
#'   \code{"top"}, \code{"left"}, \code{"right"}). Default is \code{"bottom"}.
#' @param coord_filp Logical indicating whether to flip plot coordinates using
#'   \code{coord_flip()}. Default is \code{TRUE}.
#' @param palette Integer specifying the color palette to use. Default is 3.
#' @param show_col Logical indicating whether to display color information. Default is
#'   \code{FALSE}.
#' @param cols Character vector of custom colors. If \code{NULL}, palette is used.
#'   Default is \code{NULL}.
#'
#' @return A ggplot2 object representing the stacked bar chart of cell fractions.
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD microenvironment data
#' data("sig_stad", package = "IOBR")
#' # Visualize TME cell proportions from CIBERSORT deconvolution
#' cell_bar_plot(input = sig_stad[1:20, ], id = "ID",
#'               features = colnames(sig_stad)[25:46])
cell_bar_plot <- function(input, id = "ID", title = "Cell Fraction", features = NULL, pattern = NULL, legend.position = "bottom",
                          coord_filp = TRUE, palette = 3, show_col = F, cols = NULL) {
  input <- as.data.frame(input)
  colnames(input)[which(colnames(input) == id)] <- "ID"

  if (is.null(features)) {
    if (is.null(pattern)) stop(">>>=== The 'pattern' parameter must be defined...")
    feas <- colnames(input)[str_detect(colnames(input), pattern = pattern)]
  } else {
    feas <- features
  }

  input <- input[, c("ID", feas)]

  input <- remove_names(input_df = input, variable = "colnames", patterns_to_na = patterns_to_na, patterns_space = "_")
  ##################
  if (legend.position == "top" | legend.position == "bottom") {
    legend.direction <- "horizontal"
  } else {
    legend.direction <- "vertical"
  }


  if (!is.null(cols)) {
    cols <- cols
  } else {
    if (is.null(palette)) {
      cols <- IOBR::palettes(category = "random", palette = 4, show_col = show_col, show_message = T)
    } else {
      cols <- IOBR::palettes(category = "random", palette = palette, show_col = show_col, show_message = T)
    }
  }


  if (coord_filp) {
    pp <- input %>%
      tidyr::gather(cell_type, fraction, -ID) %>%
      # plot as stacked bar chart
      ggplot(aes(x = ID, y = fraction, fill = cell_type)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      theme_light() +
      scale_fill_manual(values = cols) +
      scale_x_discrete(limits = rev(levels(input))) +
      ggtitle(paste0(title)) +
      theme(
        plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.text.x = element_text(face = "plain", angle = 0, hjust = 1, color = "black"),
        axis.text.y = element_text(face = "plain", angle = 30, hjust = 1, color = "black")
      ) +
      theme(
        legend.title = element_blank(),
        legend.position = legend.position,
        legend.direction = legend.direction,
        legend.justification = c(.5, .5),
        legend.box = "horizontal",
        legend.box.just = "top"
      )
  } else {
    pp <- input %>%
      tidyr::gather(cell_type, fraction, -ID) %>%
      # plot as stacked bar chart
      ggplot(aes(x = ID, y = fraction, fill = cell_type)) +
      geom_bar(stat = "identity") +
      # coord_flip() +
      theme_light() +
      scale_fill_manual(values = cols) +
      scale_x_discrete(limits = rev(levels(input))) +
      ggtitle(paste0(title)) +
      theme(
        plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.text.x = element_text(face = "plain", angle = 0, hjust = 1, color = "black"),
        axis.text.y = element_text(face = "plain", angle = 30, hjust = 1, color = "black")
      ) +
      theme(
        legend.title = element_blank(),
        legend.position = legend.position,
        legend.direction = legend.direction,
        legend.justification = c(.5, .5),
        legend.box = "horizontal",
        legend.box.just = "top"
      )
  }

  print(pp)
  return(pp)
}
