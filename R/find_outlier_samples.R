





#' Identification of Potential Outlier Samples in Gene Expression Data
#'
#' @description The "find_outlier_samples" function is designed to analyze gene expression data and identify potential outlier samples based on connectivity analysis. By utilizing the "WGCNA" package, this function calculates the normalized adjacency and connectivity z-scores for each sample. It further generates a connectivity plot, highlighting samples with connectivity z-scores greater than the specified y-intercept value. This function also allows for the option to plot hierarchical clustering and save the output files in a designated project folder. The returned result is a list of potential outlier samples, providing valuable insights for further analysis and data interpretation.
#'
#' @param eset A gene expression matrix data. It is the input data on which the function will operate.
#' @param yinter A numeric value representing the y-intercept for the horizontal line on the connectivity plot. It is used to identify potential outliers in the data.
#' @param project A string indicating the project name associated with the analysis. It is used to create a folder for saving the output files.
#' @param plot_hculst A logical value indicating whether to plot the hierarchical clustering of samples. If set to TRUE, the hierarchical clustering plot will be generated.
#' @param show_plot A logical value indicating whether to display the connectivity plot. If set to TRUE, the connectivity plot will be shown.
#' @param index default is null.
#'
#' @return A vector of character strings representing the names of potential outlier samples identified based on the connectivity analysis. These samples have connectivity z-scores greater than the absolute value of yinter.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # loading expression data
#' data("eset_tme_stad", package = "IOBR")
#' outs <- find_outlier_samples(eset = eset_tme_stad)
#' print(outs)
find_outlier_samples<- function(eset, yinter = -3, project = "find_outlier_eset", plot_hculst = FALSE, show_plot = TRUE, index = NULL){

  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required but not installed.")
  }
  path<- creat_folder(project)

  if(is.null(index)) index <- 1

  if(plot_hculst){
    tree.combat <- eset %>% t %>% dist %>% hclust(method = "average")
    ##############################
    pdf(paste0(path$abspath, index, "-1-clusteringplot.pdf"), width = 20, height = 10)
    plot(tree.combat, main =paste0("1-","Hierarchical Clustering Sammples"))
    dev.off()
  }
  ###############################
  ###############################
  normalized.adjacency <- (0.5 + 0.5 * bicor(eset)) ^ 2
  network.summary <- fundamentalNetworkConcepts(normalized.adjacency)
  connectivity <- network.summary$Connectivity
  connectivity.zscore <- (connectivity - mean(connectivity)) / sqrt(var(connectivity))
  connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore),
                                  Z.score = connectivity.zscore,
                                  Sample.Num = 1:length(connectivity.zscore))
  ################################

  p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Name)) +
    geom_text(size = 4, colour = "red")
  p <- p + geom_hline(aes(yintercept = yinter))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
  p<- p + design_mytheme(axis_angle = 0, axis_text_size = 12, axis_title_size = 2)
  if(show_plot) print(p)

  ggsave(p, filename = paste0(index, "-2-connectivityplot.pdf"), width = 8, height = 8, path = path$folder_name)

  names_eset_rmout <- colnames(eset)[abs(connectivity.zscore) > abs(yinter)]

  message(paste0(">>>--- When yinter = ", yinter))
  message(">>>--- Potential outliers: ")
  print(names_eset_rmout)

  return(names_eset_rmout)
}



