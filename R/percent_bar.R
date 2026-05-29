#' Create a Percent Bar Plot
#'
#' @description
#' Generates a bar plot visualizing the percentage distribution of a variable
#' grouped by another variable.
#'
#' @param input Input data frame.
#' @param x Name of the x-axis variable.
#' @param y Name of the y-axis (grouping) variable.
#' @param subset.x Optional subset of x-axis values.
#' @param color Optional color palette.
#' @param palette Optional palette type.
#' @param title Optional plot title.
#' @param axis_angle Angle for axis labels (0-90). Default is 0.
#' @param coord_flip Logical to flip coordinates. Default is FALSE.
#' @param add_Freq Logical to add frequency count. Default is TRUE.
#' @param Freq Name of frequency column.
#' @param size_freq Size of frequency labels. Default is 8.
#' @param legend.size Size of legend. Default is 0.5.
#' @param legend.size.text Size of legend text. Default is 10.
#' @param add_sum Logical to add sum to x-axis labels. Default is TRUE.
#' @param print_result Logical to print result data frame. Default is TRUE.
#' @param round.num Decimal places for proportion. Default is 2.
#'
#' @return A ggplot object.
#'
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Simulate data
#' set.seed(123)
#' sim_data <- data.frame(
#'   Subtype = sample(c("EBV", "GS", "MSI", "CIN"), 100, replace = TRUE),
#'   Lauren = sample(c("Diffuse", "Intestinal", "Mixed"), 100, replace = TRUE)
#' )
#' 
#' # Create percent bar plot
#' p <- percent_bar_plot(
#'   input = sim_data, x = "Subtype", y = "Lauren",
#'   axis_angle = 60
#' )
#' if (!is.null(p)) print(p)
percent_bar_plot <- function(input, x, y,
                             subset.x = NULL,
                             color = NULL,
                             palette = NULL,
                             title = NULL,
                             axis_angle = 0,
                             coord_flip = FALSE,
                             add_Freq = TRUE,
                             Freq = "Proportion",
                             size_freq = 8,
                             legend.size = 0.5,
                             legend.size.text = 10,
                             add_sum = TRUE,
                             print_result = TRUE,
                             round.num = 2) {
  if (is.null(input)) return(NULL)
  # Check if required columns exist
  input <- as.data.frame(input[, colnames(input) %in% c(x, y)])
  rlang::check_installed("scales")

  if (!is.null(subset.x)) {
    input <- input[input[, x] %in% subset.x, ]
  }

  if (add_Freq) {
    input$Freq <- 1
  } else {
    colnames(input)[which(colnames(input) == Freq)] <- "Freq"
  }

  df_sum <- input %>%
    dplyr::group_by(!!rlang::sym(x), !!rlang::sym(y)) %>%
    dplyr::summarise(Freq = sum(.data$Freq), .groups = "drop_last") %>%
    dplyr::group_by(!!rlang::sym(x)) %>%
    dplyr::mutate(
      Prop = round(.data$Freq / sum(.data$Freq), round.num),
      count = round(sum(.data$Freq), 0)
    )

  if (print_result) message(paste(capture.output(df_sum), collapse = "\n"))

  # Get colors
  if (is.null(color)) {
    color <- if (is.null(palette)) {
      palettes(category = "random", show_col = TRUE, show_message = TRUE)
    } else {
      palettes(category = "box", palette = palette, show_col = TRUE, show_message = TRUE)
    }
  }

  if (add_sum) {
    df_sum <- as.data.frame(df_sum)
    df_sum[, 1] <- paste0(as.character(df_sum[, 1]), "\n(", df_sum$count, ")")
  }

  hjust <- if (axis_angle == 0) 0.5 else 1

  pp <- ggplot2::ggplot(df_sum, ggplot2::aes(
    x = !!rlang::sym(x), y = .data$Prop, fill = !!rlang::sym(y)
  )) +
    ggplot2::geom_bar(stat = "identity", position = "fill", width = 0.85) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(.data$Prop, suffix = "%", accuracy = 1)),
      position = ggplot2::position_stack(0.5), size = size_freq
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::xlab(NULL) +
    design_mytheme(
      axis_text_size = 10, axis_angle = axis_angle, hjust = hjust,
      legend.size = legend.size, legend.size.text = legend.size.text
    )

  if (coord_flip) {
    pp <- pp + ggplot2::coord_flip()
  }

  if (interactive()) print(pp)
  pp
}


#' Create Pie or Donut Charts
#'
#' @description
#' Generates a pie chart or donut chart from input data.
#'
#' @param input Input dataframe.
#' @param var Variable for the chart.
#' @param var2 Secondary variable for donut chart (type = 3).
#' @param type Chart type: 1 (pie), 2 (donut), 3 (PieDonut via webr).
#' @param show_freq Logical to show frequencies. Default is FALSE.
#' @param color Optional color palette.
#' @param palette Color palette name. Default is "jama".
#' @param title Plot title. Default is NULL.
#' @param text_size Text size. Default is 10.
#' @param title_size Title size. Default is 20.
#' @param add_sum Logical to add sum to labels. Default is FALSE.
#'
#' @return A ggplot object.
#'
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Simulate data
#' set.seed(123)
#' sim_data <- data.frame(
#'   Subtype = sample(c("EBV", "GS", "MSI", "CIN"), 100, replace = TRUE)
#' )
#' 
#' # Create pie chart
#' p1 <- pie_chart(input = sim_data, var = "Subtype", palette = "jama")
#' if (!is.null(p1)) print(p1)
#' 
#' # Create donut chart
#' p2 <- pie_chart(input = sim_data, var = "Subtype", type = 2)
#' if (!is.null(p2)) print(p2)
pie_chart <- function(input, var, var2 = NULL, type = 2,
                      show_freq = FALSE, color = NULL, palette = "jama",
                      title = NULL, text_size = 10, title_size = 20,
                      add_sum = FALSE) {
  if (is.null(input)) return(NULL)
  # Check if variable exists
  input <- input[!is.na(input[, var]), , drop = FALSE]
  input <- as.data.frame(input)
  input[, var] <- as.character(input[, var])

  input2 <- input
  input <- as.data.frame(table(input[, var]))
  colnames(input)[1] <- "var"

  input <- input %>%
    dplyr::mutate(
      percent_weight = round(.data$Freq / sum(.data$Freq) * 100, 1)
    ) %>%
    dplyr::arrange(dplyr::desc(.data$Freq)) %>%
    dplyr::mutate(lab.ypos = cumsum(.data$percent_weight) - 0.5 * .data$percent_weight)

  if (interactive()) print(input)

  if (add_sum) {
    input[, 1] <- paste0(as.character(input[, 1]), "(", input$Freq, ")")
  }

  if (is.null(color)) {
    color <- if (is.null(palette)) {
      palettes(category = "random", show_col = FALSE, show_message = TRUE)
    } else {
      palettes(category = "box", palette = palette, show_col = FALSE, show_message = TRUE)
    }
  }

  pp <- switch(as.character(type),
    "1" = .pie_chart_basic(input, color, title, title_size),
    "2" = .pie_chart_donut(input, color, title, text_size, show_freq),
    "3" = .pie_chart_piedonut(input2, var, var2),
    cli::cli_abort("type must be 1, 2, or 3")
  )

  if (type %in% c(1, 2) && interactive()) print(pp)

  pp
}

#' @keywords internal
.pie_chart_basic <- function(input, color, title, title_size) {
  input <- input[order(input$Freq, decreasing = FALSE), ]

  ggplot2::ggplot(input, ggplot2::aes(x = 2, y = .data$percent_weight, fill = .data$var)) +
    ggplot2::geom_bar(stat = "identity", color = "white") +
    ggplot2::coord_polar(theta = "y", start = 0, direction = 1) +
    ggplot2::geom_text(
      ggplot2::aes(x = 2, y = .data$lab.ypos, label = .data$percent_weight),
      color = "white", size = 10
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::theme_void() +
    ggplot2::xlim(0.5, 2.5) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, color = "#666666", size = title_size)
    )
}

#' @keywords internal
.pie_chart_donut <- function(input, color, title, text_size, show_freq) {
  label_col <- if (show_freq) "Freq" else "percent_weight"
  label_fun <- function(x) if (show_freq) x else paste0(x, "%")

  ggplot2::ggplot(input, ggplot2::aes(x = 2, y = .data$percent_weight, fill = .data$var)) +
    ggplot2::geom_bar(stat = "identity", color = "white") +
    ggplot2::coord_polar(theta = "y", start = 0) +
    ggplot2::geom_text(
      ggplot2::aes(label = label_fun(!!rlang::sym(label_col))),
      color = "white", size = text_size, position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, color = "#666666", size = 20)
    )
}

#' @keywords internal
.pie_chart_piedonut <- function(input, var, var2) {
  if (is.null(var2)) {
    cli::cli_abort("var2 must be defined for type = 3")
  }
  rlang::check_installed("webr", reason = "to create PieDonut plots")

  webr::PieDonut(input, ggplot2::aes(pies = !!rlang::sym(var), donuts = !!rlang::sym(var2)),
    explode = 1, pieLabelSize = 7, donutLabelSize = 5
  )
}
