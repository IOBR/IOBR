#' Identification of TME cluster
#'
#' This function is designed to perform TME (Tumor Microenvironment) clustering analysis using the specified input dataset. The function allows for customization of various parameters to conduct clustering analysis, including the selection of features, scaling options, clustering method, and the range of cluster numbers to evaluate.
#'
#' @param input  A data frame containing the input dataset for TME clustering analysis.
#' @param features A vector specifying the features (variables) to be used for clustering. Default is NULL, which uses all columns as features.
#' @param pattern A regular expression pattern for selecting features based on column names. Default is NULL.
#' @param id A character string specifying the column in the input data frame to be used as identifiers. Default is NULL, which uses row names as identifiers.
#' @param method  A character string specifying the clustering method to use. Default is "kmeans".
#' @param min_nc An integer specifying the minimum number of clusters to evaluate. Default is 2.
#' @param max.nc An integer specifying the maximum number of clusters to evaluate. Default is 6.
#' @param scale A logical value indicating whether to scale the selected features. Default is TRUE.
#'
#' @return A data frame with clusters appended as a new column.
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data("tcga_stad_sig", package = "IOBR")
#' res <- tme_cluster(input = tcga_stad_sig, features = NULL, pattern = "xCell", id = "ID", method = "kmeans", min_nc = 2, max.nc = 6)
#' sig_heatmap(input = res, features = colnames(res)[3:ncol(res)], group = "cluster")
tme_cluster <- function(input, features = NULL, pattern = NULL, id = NULL, scale = TRUE, method = "kmeans", min_nc = 2, max.nc = 6) {
  input <- as.data.frame(input)

  if (is.null(id)) {
    id <- rownames(input)
  } else {
    id <- as.character(input[, id])
  }

  if (is.null(features)) {
    if (!is.null(pattern)) {
      features <- colnames(input)[str_detect(colnames(input), pattern = pattern)]
    } else {
      features <- colnames(input)
    }
    input <- input[, colnames(input) %in% features]
  } else {
    features <- features[features %in% colnames(input)]
    input <- input[, colnames(input) %in% features]
  }

  if (scale) {
    input[, features] <- scale(input[, features], scale = TRUE, center = TRUE)
    features <- feature_manipulation(data = input, feature = features, print_result = FALSE)
    input <- input[, colnames(input) %in% features]
  }
  # print(input)
  res <- NbClust::NbClust(
    data = input,
    diss = NULL,
    distance = "euclidean",
    min.nc = min_nc,
    max.nc = max.nc,
    method = method,
    index = "kl",
    alphaBeale = 0.1
  )
  nc <- res$Best.nc

  message(print(">>>== Best number of TME clusters is: "))
  print(nc)

  out <- res$Best.partition

  message(print(">>>== Cluster of samples: "))

  out <- data.frame("ID" = id, "cluster" = out)
  out$cluster <- paste0("TME", out$cluster)
  print(summary(as.factor(out$cluster)))

  input$ID <- id
  res <- merge(out, input, by = "ID")
  # print(out)
  return(res)
}
