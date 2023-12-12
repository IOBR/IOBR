



#' Principal Component Analysis (PCA) Visualization Function
#' @description The iobr_pca function performs Principal Component Analysis (PCA), which reduces the dimensionality of data while maintaining most of the original variance, and visualizes the PCA results on a scatter plot.
#'
#' @param data The input data for PCA. It should be a matrix or a data frame.
#' @param is.matrix Specifies whether the input data is a matrix. Default is TRUE.
#' @param scale Specifies whether to scale the input data. Default is TRUE.
#' @param is.log Specifies whether to log transform the input data. Default is FALSE.
#' @param pdata Additional data associated with the principal components. It should be a data frame. Default is NULL.
#' @param id_pdata The column name in 'pdata' that represents the ID for matching with 'data'. Default is "ID".
#' @param group  The column name in 'pdata' that represents groups/categories to color the points. Default is NULL.
#' @param cols The color scheme to be used for group categories. Default is "normal".
#' @param palette The color palette to be used for group categories. Default is "jama".
#' @param repel Specifies whether to repel the data points to avoid overlap. Default is FALSE.
#' @param ncp The number of dimensions to keep in the PCA. Default is 5.
#' @param axes  The dimensions/axes to be plotted. Default is c(1, 2).
#' @param addEllipses Specifies whether to add concentration ellipses to the plot. Default is TRUE.
#' @param geom.ind The type of geometric representation for the points in the PCA plot. Default is "point".
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data("eset_stad", package = "IOBR")
#' eset <- count2tpm(eset_stad)
#' iobr_pca(eset, is.matrix = TRUE, scale = TRUE, is.log = TRUE, pdata = stad_group, id_pdata = "ID", group = "subtype")
#'
iobr_pca <- function(data, is.matrix = TRUE, scale = TRUE, is.log = FALSE, pdata, id_pdata = "ID", group = NULL,
                     geom.ind = "point", cols = "normal", palette = "jama", repel = FALSE, ncp = 5, axes = c(1, 2), addEllipses = TRUE){

  if(is.log) data <- log2eset(data+1)
  feas <- feature_manipulation(data = data, feature = rownames(data), is_matrix = TRUE)
  data <-data[rownames(data)%in%feas, ]
  #######################################
  if(is.matrix) data <- t(data)
  if(scale) data <- scale(data)

  res.pca <- FactoMineR::PCA(data, ncp = ncp, graph = FALSE)

  pdata <- as.data.frame(pdata)
  colnames(pdata)[which(colnames(pdata)==id_pdata)] <- "id"

  pdata <- pdata[pdata$id%in%rownames(data), ]
  pdata <- pdata[match(rownames(data), pdata$id), ]

  message(print(table(pdata[, group])))
  ##########################################

  cols <- get_cols(cols = cols, palette = palette, show_col = FALSE, seed = 123)

  # print(cols)
  # print(pdata[, group])
  cols <- cols[1:length(unique(pdata[, group]))]

  print(paste0(">>-- colors for PCA: "))
  message(paste0(">>== ", cols))
  #########################################
  p <- factoextra:: fviz_pca_ind(res.pca,
                                 axes = axes,
                                 geom.ind     = geom.ind, # show points only (nbut not "text")
                                 col.ind      = as.character(pdata[, group]), # color by groups
                                 palette      = cols,
                                 repel        = repel,
                                 addEllipses  = addEllipses, # Concentration ellipses
                                 legend.title = group)
  return(p)

}
