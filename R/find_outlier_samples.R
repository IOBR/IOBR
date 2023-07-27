





#' Title
#'
#' @param eset expression set
#' @param yinter between -3 to -2
#' @param project default is NSCLC
#' @param plot_hculst
#'
#' @return
#' @export
#'
#' @examples
find_outlier_samples<- function(eset, yinter = -3, project = "NSCLC", plot_hculst = FALSE){

  library(WGCNA)
  path<- creat_folder(project)


  if(plot_hculst){
    tree.combat <- eset %>% t %>% dist %>% hclust(method = "average")
    ##############################
    pdf(paste0(path$abspath, "1-clusteringplot.pdf"), width = 20, height = 10)
    plot(tree.combat, main =paste0("1-","Hierarchical Clustering Sammples"))
    dev.off()
  }
  ###############################

  #' 使用WGCNA查看样本离异度
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
  p<-p+design_mytheme()

  ggsave(p, filename = paste0("2-connectivityplot.pdf"), width = 10, height = 10, path = path$folder_name)

  names_eset_rmout <- colnames(eset)[abs(connectivity.zscore) > abs(yinter)]

  message(paste0(">>>--- when yinter = ", yinter))
  message(">>>--- Potential outliers: ")
  print(names_eset_rmout)

  return(names_eset_rmout)
}


# find_outlier_samples(eset = eset, project = "melanoma")

