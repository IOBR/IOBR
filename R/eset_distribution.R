#' Title eset_distribution
#'
#' @description eset_distribution generates boxplots and density plots to analyze the distribution of expression values in an expression set (eset).
#' @param eset Expression set (matrix or data frame) containing gene expression data.
#' @param quantile Number of quantiles used for sampling columns in the eset matrix.
#' @param log Logical value indicating whether to perform log2 transformation of the expression values.
#' @param project Optional string specifying the name of the project. If not provided, a default name ('ESET') is used.
#'
#' @return This function does not return a value but saves boxplots and density plots as PNG files in the specified or default directory.
#' @export
#'
#' @examples
#' data("eset_stad", package = "IOBR")
#' eset <- anno_eset(eset = eset_stad, annotation = anno_rnaseq)
#' eset_distribution(eset)
eset_distribution <- function(eset, quantile = 3, log = TRUE, project = NULL) {
  feas <- feature_manipulation(data = eset, feature = rownames(eset), is_matrix = TRUE)
  eset <- eset[rownames(eset) %in% feas, ]

  index <- colnames(eset)[sample(1:ncol(eset), round(dim(eset)[2] / quantile, 0))]
  eset1 <- eset[, index]
  eset1 <- log2eset(eset1)

  eset1 <- eset1 %>%
    t() %>%
    data.frame(Sample.Name = colnames(eset1))
  eset.melt <- reshape2::melt(eset1, id = c("Sample.Name"))
  head(eset.melt)
  colnames(eset.melt)[2:3] <- c("Symbol", "Intensity")
  ######################################
  if (is.null(project)) {
    path <- creat_folder("result")
    project <- "ESET"
  } else {
    path <- creat_folder(project)
  }

  p <- ggplot(eset.melt, aes(x = Sample.Name, y = Intensity)) +
    geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(size = 23))
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0))
  p <- p + ggtitle(
    "Normalized signal intensity",
    paste0(
      "Patient No. is:  ", dim(eset)[2], "; ",
      "Number of features is: ", dim(eset)[1], ";  ",
      "Maximum is: ", round(max(eset), 2)
    )
  ) +
    ylab("Intensity") + xlab("Sample")
  p <- p + theme(
    legend.position = "none", panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
    theme(plot.subtitle = element_text(size = 18, hjust = 0.1, face = "italic", color = "black"))
  ###################################
  ggsave(
    filename = paste0("1-", project, "-boxplot.png"), plot = p,
    width = 15, height = 8, path = path$folder_name
  )
  ####################################

  p <- ggplot(eset.melt, aes(Intensity, group = Sample.Name)) +
    geom_density() +
    theme_bw()
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + ggtitle("Histogram of Log2 Expression") + ylab("Density") + xlab("Log2 Expression")
  ggsave(
    filename = paste0("2-", project, "-Densityplot.png"), plot = p,
    width = 9, height = 6, path = path$folder_name
  )
}
