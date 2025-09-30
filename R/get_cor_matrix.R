#' Calculate and Visualize Correlation Matrix between Two Variable Sets
#' @description The "get_cor_matrix" function calculates and visualizes the correlation matrix between two sets of variables in a dataset. It provides flexibility in defining correlation methods, handling missing values, and incorporating additional data. The function supports various correlation methods, such as Pearson correlation, and displays the correlation result in a customizable plot. The plot includes color-coded tiles representing the correlation value between each pair of variables. Additionally, the function offers options to save the plot and data for further analysis.
#' @param data The input dataset for correlation analysis. Genes or signatures in columns. If data is a matrix, please set `is.matrix` = TRUE, and data will be t()
#' @param feas1  A vector of variable names representing the first set of variables for correlation analysis.
#' @param feas2 A vector of variable names representing the second set of variables for correlation analysis.
#' @param method The method used to calculate correlation. Default is "pearson".
#' @param path The path to save the correlation plot. If not specified, the plot will be saved in a folder named "index-cor_matrix_plot".
#' @param index A numeric value used for naming the output plot file. Default is 1.
#' @param fig.type The file format of the output plot. Default is "pdf".
#' @param width The width of the output plot.
#' @param height The height of the output plot.
#' @param project The name of the project or analysis. This will be used as the plot title. Default is NULL.
#' @param is.matrix A logical value indicating whether the input data is a matrix. Default is FALSE.
#' @param scale A logical value indicating whether to scale the data before correlation analysis. Default is TRUE.
#' @param font.size The font size of the axis labels and tick marks in the plot. Default is 15.
#' @param fill_by_cor A logical value indicating whether to fill the tiles in the plot according to the correlation values. Default is FALSE.
#' @param round.num The number of decimal places to round the correlation values. Default is 1.
#' @param font.size.star The font size of the significance stars in the plot. Default is 8.
#'
#' @return A ggplot object displaying the correlation matrix.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Assuming 'data' is a data frame or matrix with genes or signatures as columns
#' set.seed(123)
#' data <- as.data.frame(matrix(rnorm(1000), nrow = 100, ncol = 10))
#' colnames(data) <- paste("Gene", 1:10, sep = "_")
#'
#' # Define two sets of variables for correlation analysis
#' feas1 <- c("Gene_1", "Gene_2", "Gene_3")
#' feas2 <- c("Gene_4", "Gene_5", "Gene_6")
#'
#' # Generate and visualize the correlation matrix
#' cor_plot <- get_cor_matrix(
#'   data = data, feas1 = feas1, feas2 = feas2, scale = TRUE,
#'   method = "spearman", project = "Example Correlation Matrix"
#' )
#' print(cor_plot)
get_cor_matrix <- function(data, feas1, feas2, method = "pearson", path = NULL, index = 1, fig.type = "pdf", width = NULL, height = NULL,
                           project = NULL, is.matrix = FALSE, scale = T, font.size.star = 8, font.size = 15, fill_by_cor = FALSE, round.num = 1) {
  if (is.null(path)) {
    path <- creat_folder(paste0(index, "-cor_matrix_plot"))
  } else {
    path <- creat_folder(path)
  }

  if (is.matrix) data <- as.data.frame(t(data))

  feas1 <- feas1[feas1 %in% colnames(data)]
  feas2 <- feas2[feas2 %in% colnames(data)]

  if (scale) data <- scale(data[, unique(c(feas1, feas2))])
  #############################
  result <- psych::corr.test(data[, feas1], data[, feas2], method = method)
  #############################
  heat <- cbind(reshape2::melt(result$r), reshape2::melt(result$p))[, c(1, 2, 3, 6)]
  colnames(heat) <- c("ID1", "ID2", "cor", "pvalue")
  #############################

  if (fill_by_cor) {
    heat$stars <- round(heat$cor, round.num)
  } else {
    heat$stars <- cut(heat$pvalue,
      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.5, Inf),
      label = c("***", "**", "*", "+", "")
    )
  }

  heat$ID2 <- gsub(heat$ID2, pattern = "\\_", replacement = " ")
  #######################################
  # Plot everything
  p <- ggplot(aes(x = ID1, y = ID2, fill = cor), data = heat)
  ########################################
  fig12 <- p + geom_tile() + scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +
    # geom_text(aes(label=stars, color=value), size=8) +
    # scale_colour_gradient(low="grey30", high="white", guide="none") +
    geom_text(aes(label = stars), color = "black", size = font.size.star) +
    scale_y_discrete(name = "signature", labels = function(y) str_wrap(y, width = 40)) +
    labs(y = NULL, x = NULL, fill = "Coefficient") +
    # geom_vline(xintercept=1.5, size=1.5, color="grey50") +
    theme_bw() + theme(axis.text.x = element_text(angle = -45, hjust = 0, size = font.size)) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = font.size)) +
    theme(axis.title = element_text(angle = 0, hjust = 0.5, size = font.size + 4))

  if (is.null(project)) {
    project <- "-"
  } else {
    fig12 <- fig12 + ggtitle(label = project)
  }

  #########################################
  if (is.null(width)) width <- length(feas1) * 0.55 + 6.5
  if (is.null(height)) height <- length(feas2) * 0.35 + 4.5

  ggsave(fig12,
    filename = paste0(index, "-", project, "-cor_plot.", fig.type),
    width = width, height = height, path = path$folder_name
  )
  print(fig12)
  #######################################################
  return(fig12)
}
