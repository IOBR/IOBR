#' Create a Percent Bar Plot
#'
#' This function generates a bar plot that visualizes the percentage distribution of a variable grouped by another variable.
#' It offers options to customize the plot's colors, title, axis labels, and more.
#'
#' @param input Input data frame.
#' @param x Name of the x-axis variable.
#' @param y Name of the y-axis variable.
#' @param color Optional color palette for the bars.
#' @param title Optional title for the plot.
#' @param palette Optional palette type for color selection.
#' @param axis_angle Optional angle for the axis labels. range from [0,90]
#' @param add_Freq default is true for data frame
#' @param Freq Optional name for the frequency column.
#' @param size_freq Size of the frequency labels.
#' @param add_sum Boolean indicating whether to add the sum to the x-axis labels.
#' @param subset.x Optional subset of x-axis values.
#' @param legend.size Size of the legend.
#' @param legend.size.text Size of the legend text.
#' @param print_result Boolean indicating whether to print the result data frame.
#' @param round.num Number of decimal places to round the proportion column.
#' @param coord_flip Boolean indicating whether to flip the x and y axes.
#'
#' @return A ggplot object representing the percentage bar plot.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("sig_stad", package = "IOBR")
#' table(sig_stad$Subtype, sig_stad$Lauren)
#' percent_bar_plot(input = sig_stad, x = "Subtype", y = "Lauren", axis_angle = 60)
percent_bar_plot <- function(input, x, y,
                             subset.x = NULL,
                             color = NULL,
                             palette = NULL,
                             title = NULL,
                             axis_angle = 0,
                             coord_flip = FALSE,
                             add_Freq = TRUE,
                             Freq = NULL,
                             size_freq = 8,
                             legend.size = 0.5,
                             legend.size.text = 10,
                             add_sum = T,
                             print_result = T,
                             round.num = 2) {
  input <- as.data.frame(input[, colnames(input) %in% c(x, y)])

  if (!is.null(subset.x)) {
    input <- input[input[, x] %in% subset.x, ]
  }

  if (add_Freq) {
    input$Freq <- 1
  } else {
    # input$Freq<-1
    colnames(input)[which(colnames(input) == Freq)] <- "Freq"
  }
  # bpercent<-as.data.frame(table(input[,x],input[,y]))
  # colnames(bpercent)<-c(x,y,"Frequency")
  # ce =plyr::ddply(bpercent, y, transform, percent_weight = Frequency / sum(Frequency))
  # ce<-ce[order(ce[,y],decreasing = F),]
  # print(ce)


  df_sum <- input %>%
    dplyr::group_by(!!sym(x), !!sym(y)) %>%
    dplyr::summarise(Freq = sum(Freq)) %>%
    dplyr::group_by(!!sym(x)) %>%
    dplyr::mutate(Prop = round(Freq / sum(Freq), round.num)) %>%
    dplyr::mutate(count = round(sum(Freq), 0))

  if (print_result) print(df_sum)

  # print(as.data.frame(df_sum))

  if (!is.null(color)) {
    color <- color
  } else {
    if (is.null(palette)) {
      color <- IOBR::palettes(category = "random", show_col = T, show_message = T)
    } else {
      color <- IOBR::palettes(category = "box", palette = palette, show_col = T, show_message = T)
    }
  }


  if (add_sum) {
    df_sum <- as.data.frame(df_sum)
    df_sum[, 1] <- paste0(as.character(df_sum[, 1]), "(", df_sum$count, ")")
  }

  # https://github.com/tidyverse/ggplot2/issues/3369
  c <- ggplot(df_sum, aes(x = !!sym(x), y = Prop, fill = !!sym(y))) +
    geom_bar(stat = "identity", position = "fill", width = 0.85) +
    geom_text(aes(label = scales::percent(Prop, suffix = "%", accuracy = 1)), position = position_stack(.5), size = size_freq) + #
    ggtitle(title) +
    scale_fill_manual(values = color) +
    xlab(NULL)

  if (axis_angle == 0) {
    mytheme <- design_mytheme(axis_text_size = 20, axis_angle = 0, hjust = 0.5, legend.size = legend.size, legend.size.text = legend.size.text)
  } else {
    mytheme <- design_mytheme(axis_text_size = 20, axis_angle = axis_angle, hjust = 1, legend.size = legend.size, legend.size.text = legend.size.text)
    message(">>>=== When the coordinates are inverted, axis_angle can't fulfil its function")
  }

  pp <- c + mytheme
  if (coord_flip) {
    pp <- pp + coord_flip()
  }
  print(pp)
  return(pp)
}






#' pie_chart
#'
#' @description
#' This function generates a pie chart or a donut chart from the input data. It allows customization of various visual aspects such as colors, labels, and title.
#'
#' @param input The input dataframe.
#' @param var The variable on which the pie chart or donut chart will be based.
#' @param color Optional. The color palette for the chart. If not provided, a default color palette will be used.
#' @param palette Optional. The color palette to be used if color is not specified. Default is "jama".
#' @param title Optional. The title of the chart. Default is NULL.
#' @param text_size  The size of the text on the chart. Default is 10.
#' @param title_size The size of the title text on the chart. Default is 20.
#' @param var2 Optional. A secondary variable for creating a donut chart. Default is NULL.
#' @param type The type of chart to be generated. 1 for pie chart, 2 for donut chart, 3 for donut chart based on webr package.
#' @param show_freq Boolean indicating whether to show frequencies on the chart. Default is FALSE.
#' @param add_sum Boolean indicating whether to add the sum of frequencies to the chart. Default is FALSE.
#'
#' @return Generate Pie or Donut Charts
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' pie_chart(input = sig_stad, var = "Subtype", palette = "jama")
#' pie_chart(input = sig_stad, var = "Subtype", palette = "nrc")
pie_chart <- function(input, var, var2 = NULL, type = 2, show_freq = FALSE, color = NULL, palette = "jama", title = NULL, text_size = 10, title_size = 20, add_sum = FALSE) {
  input <- input[!is.na(input[, var]), ]
  input <- as.data.frame(input)
  input[, var] <- as.character(input[, var])

  input2 <- input
  input <- as.data.frame(table(input[, var]))

  colnames(input)[1] <- "var"

  input <- input %>%
    mutate(percent_weight = round(Freq / sum(Freq) * 100, 1)) %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    dplyr::mutate(lab.ypos = cumsum(percent_weight) - 0.5 * percent_weight)
  print(input)


  if (add_sum) {
    input <- as.data.frame(input)
    input[, 1] <- paste0(as.character(input[, 1]), "(", input$Freq, ")")
  }


  # if(!is.null(color)) color<-IOBR::palettes(category = "random",show_col = F,show_message = F)
  # if(!is.null(palette)) color<-IOBR::palettes(category = "box",palette = palette, show_col = F,show_message = F)

  if (!is.null(color)) {
    color <- color
  } else {
    if (is.null(palette)) {
      color <- IOBR::palettes(category = "random", show_col = T, show_message = T)
    } else {
      color <- IOBR::palettes(category = "box", palette = palette, show_col = T, show_message = T)
    }
  }


  if (type == 1) {
    input <- input[order(input$Freq, decreasing = F), ]
    pp <- ggplot(input, aes(x = 2, y = percent_weight, fill = var)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar(theta = "y", start = 0, direction = 1) +
      geom_text(aes(x = 2, y = lab.ypos, label = percent_weight), color = "white", size = 10) +
      labs(x = NULL, y = NULL, fill = NULL) +
      ggtitle(paste0(title)) +
      scale_fill_manual(values = color) +
      theme_void() +
      xlim(0.5, 2.5) +
      theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666", size = title_size)
      )
  }

  if (type == 2) {
    if (show_freq) {
      pp <- ggplot(input, aes(x = 2, y = percent_weight, fill = var)) +
        geom_bar(stat = "identity", color = "white") +
        coord_polar(theta = "y", start = 0) +
        geom_text(aes(label = paste0(Freq)),
          color = "white",
          size = text_size, position = position_stack(vjust = 0.5)
        ) +
        labs(x = NULL, y = NULL, fill = NULL) +
        ggtitle(paste0(title)) +
        scale_fill_manual(values = color) +
        theme_classic() +
        theme(
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, color = "#666666", size = title_size)
        )
    } else {
      pp <- ggplot(input, aes(x = 2, y = percent_weight, fill = var)) +
        geom_bar(stat = "identity", color = "white") +
        coord_polar(theta = "y", start = 0) +
        geom_text(aes(label = paste0(percent_weight, "%")),
          color = "white",
          size = text_size, position = position_stack(vjust = 0.5)
        ) +
        labs(x = NULL, y = NULL, fill = NULL) +
        ggtitle(paste0(title)) +
        scale_fill_manual(values = color) +
        theme_classic() +
        theme(
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, color = "#666666", size = title_size)
        )
    }
  }

  if (type == 3) {
    # https://rpubs.com/cardiomoon/398623
    if (is.null(var2)) stop("var2 must be defined!")
    pp <- webr::PieDonut(input2, aes(pies = !!sym(var), donuts = !!sym(var2)),
      explode = 1, pieLabelSize = 7,
      donutLabelSize = 5
    )
  }

  if (type == 1 | type == 2) print(pp)

  return(pp)
}
