#' Identify Outlier Samples in Gene Expression Data
#'
#' Analyzes gene expression data to identify potential outlier samples using connectivity analysis via the WGCNA package. Calculates normalized adjacency and connectivity z-scores for each sample, generates connectivity plots, and optionally performs hierarchical clustering. Outlier samples are those with connectivity z-scores exceeding the specified threshold.
#'
#' @param eset Numeric matrix. Gene expression data with genes as rows and samples as columns.
#' @param yinter Numeric. Y-intercept threshold for identifying outliers in the connectivity plot. Default is -3.
#' @param project Character. Project name for output folder creation. Default is "find_outlier_eset".
#' @param plot_hculst Logical. Whether to plot hierarchical clustering of samples. Default is FALSE.
#' @param show_plot Logical. Whether to display the connectivity plot. Default is TRUE.
#' @param index Integer or NULL. Index for output file naming. Default is NULL.
#'
#' @return Character vector of sample names identified as potential outliers (connectivity z-score > |yinter|).
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Load expression data
#' data("eset_tme_stad", package = "IOBR")
#' outs <- find_outlier_samples(eset = eset_tme_stad)
#' print(outs)
find_outlier_samples <- function(eset, yinter = -3, project = "find_outlier_eset", plot_hculst = FALSE, show_plot = TRUE, index = NULL) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required but not installed.")
  }
  path <- creat_folder(project)

  if (is.null(index)) index <- 1

  if (plot_hculst) {
    tree.combat <- eset %>%
      t() %>%
      dist() %>%
      hclust(method = "average")
    ##############################
    pdf(paste0(path$abspath, index, "-1-clusteringplot.pdf"), width = 20, height = 10)
    plot(tree.combat, main = paste0("1-", "Hierarchical Clustering Sammples"))
    dev.off()
  }
  ###############################
  ###############################
  normalized.adjacency <- (0.5 + 0.5 * bicor(eset))^2
  network.summary <- fundamentalNetworkConcepts(normalized.adjacency)
  connectivity <- network.summary$Connectivity
  connectivity.zscore <- (connectivity - mean(connectivity)) / sqrt(var(connectivity))
  connectivity.plot <- data.frame(
    Sample.Name = names(connectivity.zscore),
    Z.score = connectivity.zscore,
    Sample.Num = 1:length(connectivity.zscore)
  )
  ################################

  p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Name)) +
    geom_text(size = 4, colour = "red")
  p <- p + geom_hline(aes(yintercept = yinter))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
  p <- p + design_mytheme(axis_angle = 0, axis_text_size = 12, axis_title_size = 2)
  if (show_plot) print(p)

  ggsave(p, filename = paste0(index, "-2-connectivityplot.pdf"), width = 8, height = 8, path = path$folder_name)

  names_eset_rmout <- colnames(eset)[abs(connectivity.zscore) > abs(yinter)]

  message(paste0(">>>--- When yinter = ", yinter))
  message(">>>--- Potential outliers: ")
  print(names_eset_rmout)

  return(names_eset_rmout)
}
