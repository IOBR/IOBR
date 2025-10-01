#' Design Custom Theme for ggplot2 Plots
#'
#' @description This function creates a customized ggplot2 theme based on user-specified parameters for plot elements such as title size, axis sizes, legend settings, and theme style. It supports various base themes and allows fine-tuning of visual aspects for publication-quality plots.
#' @param theme The base theme for the plot. Options include "light", "bw", "classic", and "classic2". Default is "light".
#' @param plot_title_size The relative size of the plot title. Default is 2.
#' @param axis_title_size The relative size of the axis titles. Default is 2.
#' @param axis_text_size The size of the axis tick mark labels. Default is 12.
#' @param axis_angle The angle of rotation for the x-axis tick mark labels. Default is 60.
#' @param hjust Horizontal justification for x-axis text. Default is 1.
#' @param legend.position The position of the legend in the plot. Options include "none", "left", "right", "bottom", and "top". Default is "bottom".
#' @param legend.direction The direction of the legend items. Options include "horizontal" and "vertical". Default is "horizontal".
#' @param legend.size The size of the legend key. Default is 0.25.
#' @param legend.key.height The height of the legend key in cm. Default is 0.5.
#' @param legend.key.width The width of the legend key in cm. Default is 0.5.
#' @param legend.size.text The size of the legend text labels. Default is 10.
#' @param legend.box The orientation of the legend box. Options include "horizontal" and "vertical". Default is "horizontal".
#'
#' @author Dongqiang Zeng
#' @return A ggplot2 theme object that can be added to a ggplot.
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' mytheme <- design_mytheme(theme = "bw", plot_title_size = 1.5, axis_text_size = 14)
#' p + mytheme + ggtitle("Example Plot")
design_mytheme <- function(theme = "light",
                           plot_title_size = 2,
                           axis_title_size = 2,
                           axis_text_size = 12,
                           axis_angle = 60,
                           hjust = 1,
                           legend.position = NULL,
                           legend.direction = NULL,
                           legend.size = 0.25,
                           legend.key.height = 0.5,
                           legend.key.width = 0.5,
                           legend.size.text = 10,
                           legend.box = "horizontal") {
  if (theme == "light") {
    mytheme <- ggplot2::theme_light()
  } else if (theme == "bw") {
    mytheme <- ggplot2::theme_bw()
  } else if (theme == "classic") {
    mytheme <- ggplot2::theme_classic()
  } else if (theme == "classic2") {
    mytheme <- ggpubr::theme_classic2()
  }
  message(paste0(">>>>Options for `theme`: light, bw, classic and classic2"))

  mytheme <- mytheme + ggplot2::theme(
    plot.title = element_text(size = rel(plot_title_size), hjust = 0.5),
    axis.title = element_text(size = rel(axis_title_size)),
    axis.text = element_text(size = rel(4)),
    axis.text.x = element_text(
      face = "plain", size = axis_text_size,
      angle = axis_angle, hjust = hjust, color = "black"
    ), # family="Times New Roman"
    axis.text.y = element_text(face = "plain", size = axis_text_size, color = "black"), # family="Times New Roman"
    # panel.grid.major=element_line(color="white"),
    # panel.grid.minor=element_line(color="white"),
    # panel.border=element_rect(color="white"),
    axis.line = element_line(color = "grey", size = 0.2)
  )

  if (is.null(legend.position)) {
    legend.position <- "bottom"
  } else {
    message(paste0(">>>>Options for 'legend.position' : none, left, right, bottom, top"))
  }

  if (is.null(legend.direction)) {
    legend.direction <- "horizontal"
  } else {
    message(paste0(">>>>Options for 'legend.direction' : horizontal, vertical "))
  }


  mytheme <- mytheme +
    theme(
      # legend.key.size=unit(legend.size,"cm"),
      legend.key.height = unit(legend.key.height, "cm"),
      legend.key.width = unit(legend.key.width, "cm"),
      # legend.title=element_blank(),
      legend.position = legend.position, # "none","left","right","bottom","top",or #c(0.5,1)
      legend.direction = legend.direction, # "vertical"
      legend.justification = c(.5, .5), # "center" or two-element numeric vector
      legend.box = legend.box, # "horizontal",
      legend.box.just = "top",
      legend.text = element_text(colour = "black", size = legend.size.text, face = "plain") # ("plain", "italic", "bold", "bold.italic")
    )

  return(mytheme)
}
