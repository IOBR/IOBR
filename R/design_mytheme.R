#' Design Custom Theme for ggplot2 Plots
#'
#' @description
#' Creates a customized ggplot2 theme based on user-specified parameters for
#' plot elements such as title size, axis sizes, legend settings, and theme style.
#' Supports various base themes and allows fine-tuning of visual aspects.
#'
#' @param theme Base theme: "light", "bw", "classic", "classic2".
#'   Default is "light".
#' @param plot_title_size Relative size of plot title. Default is 2.
#' @param axis_title_size Relative size of axis titles. Default is 2.
#' @param axis_text_size Size of axis tick labels. Default is 12.
#' @param axis_angle Angle of x-axis tick labels. Default is 60.
#' @param hjust Horizontal justification for x-axis text. Default is 1.
#' @param legend.position Legend position: "none", "left", "right", "bottom",
#'   "top". Default is "bottom".
#' @param legend.direction Direction of legend items: "horizontal" or "vertical".
#'   Default is "horizontal".
#' @param legend.size Size of legend key. Default is 0.25.
#' @param legend.key.height Height of legend key in cm. Default is 0.5.
#' @param legend.key.width Width of legend key in cm. Default is 0.5.
#' @param legend.size.text Size of legend text labels. Default is 10.
#' @param legend.box Orientation of legend box: "horizontal" or "vertical".
#'   Default is "horizontal".
#'
#' @return A ggplot2 theme object.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point()
#' mytheme <- design_mytheme(theme = "bw", plot_title_size = 1.5, axis_text_size = 14)
#' p + mytheme + ggtitle("Example Plot")
design_mytheme <- function(theme = c("light", "bw", "classic", "classic2"),
                           plot_title_size = 2,
                           axis_title_size = 2,
                           axis_text_size = 12,
                           axis_angle = 60,
                           hjust = 1,
                           legend.position = "bottom",
                           legend.direction = "horizontal",
                           legend.size = 0.25,
                           legend.key.height = 0.5,
                           legend.key.width = 0.5,
                           legend.size.text = 10,
                           legend.box = "horizontal") {
  theme <- rlang::arg_match(theme)

  # Base theme
  base_theme <- switch(theme,
    light = ggplot2::theme_light(),
    bw = ggplot2::theme_bw(),
    classic = ggplot2::theme_classic(),
    classic2 = ggplot2::theme_classic() + ggplot2::theme(text = ggplot2::element_text(size = 14))
  )

  # Customize theme
  mytheme <- base_theme + ggplot2::theme(
    plot.title = ggplot2::element_text(size = ggplot2::rel(plot_title_size), hjust = 0.5),
    axis.title = ggplot2::element_text(size = ggplot2::rel(axis_title_size)),
    axis.text.x = ggplot2::element_text(
      face = "plain", size = axis_text_size,
      angle = axis_angle, hjust = hjust, color = "black"
    ),
    axis.text.y = ggplot2::element_text(
      face = "plain", size = axis_text_size, color = "black"
    ),
    axis.line = ggplot2::element_line(color = "grey", linewidth = 0.2)
  )

  # Legend settings
  mytheme <- mytheme + ggplot2::theme(
    legend.key.height = ggplot2::unit(legend.key.height, "cm"),
    legend.key.width = ggplot2::unit(legend.key.width, "cm"),
    legend.position = legend.position,
    legend.direction = legend.direction,
    legend.justification = c(0.5, 0.5),
    legend.box = legend.box,
    legend.box.just = "top",
    legend.text = ggplot2::element_text(
      colour = "black", size = legend.size.text, face = "plain"
    )
  )

  mytheme
}
