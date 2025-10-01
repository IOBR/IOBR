#' Principal Component Analysis (PCA) Visualization
#'
#' This function performs Principal Component Analysis (PCA) on gene expression data,
#' reduces dimensionality while preserving variance, and generates a scatter plot visualization.
#'
#' @param data Input data for PCA: matrix or data frame.
#' @param is.matrix Logical indicating if input is a matrix. Default is TRUE.
#' @param scale Logical indicating whether to scale the data. Default is TRUE.
#' @param is.log Logical indicating whether to log-transform the data. Default is FALSE.
#' @param pdata Data frame with sample IDs and grouping information.
#' @param id_pdata Column name in `pdata` for sample IDs. Default is "ID".
#' @param group Column name in `pdata` for grouping variable. Default is NULL.
#' @param cols Color scheme for groups. Default is "normal".
#' @param palette Color palette for groups. Default is "jama".
#' @param repel Logical indicating whether to repel overlapping points. Default is FALSE.
#' @param ncp Number of principal components to retain. Default is 5.
#' @param axes Principal components to plot (e.g., c(1, 2)). Default is c(1, 2).
#' @param addEllipses Logical indicating whether to add concentration ellipses. Default is TRUE.
#' @param geom.ind Type of geometric representation for points. Default is "point".
#'
#' @return A ggplot object of the PCA plot.
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data("eset_stad", package = "IOBR")
#' eset <- count2tpm(eset_stad)
#' iobr_pca(eset, is.matrix = TRUE, scale = TRUE, is.log = TRUE, pdata = stad_group, id_pdata = "ID", group = "subtype")
#'
iobr_pca <- function(data, is.matrix = TRUE, scale = TRUE, is.log = FALSE, pdata, id_pdata = "ID", group = NULL,
                     geom.ind = "point", cols = "normal", palette = "jama", repel = FALSE, ncp = 5, axes = c(1, 2), addEllipses = TRUE) {
  if (is.log) data <- log2eset(data + 1)
  feas <- feature_manipulation(data = data, feature = rownames(data), is_matrix = TRUE)
  data <- data[rownames(data) %in% feas, ]
  #######################################
  if (is.matrix) data <- t(data)
  if (scale) data <- scale(data)

  res.pca <- FactoMineR::PCA(data, ncp = ncp, graph = FALSE)

  pdata <- as.data.frame(pdata)
  colnames(pdata)[which(colnames(pdata) == id_pdata)] <- "id"

  pdata <- pdata[pdata$id %in% rownames(data), ]
  pdata <- pdata[match(rownames(data), pdata$id), ]


  message(print(table(pdata[, group])))
  ##########################################

  cols <- get_cols(cols = cols, palette = palette, show_col = FALSE, seed = 123)

  # print(cols)
  # print(pdata[, group])
  cols <- cols[1:length(unique(pdata[[group]]))]

  print(paste0(">>== colors for group: "))
  message(paste0(">>== ", cols))
  #########################################
  p <- factoextra::fviz_pca_ind(res.pca,
    axes = axes,
    geom.ind = geom.ind, # show points only (nbut not "text")
    col.ind = as.character(pdata[, group]), # color by groups
    palette = cols,
    repel = repel,
    addEllipses = addEllipses, # Concentration ellipses
    legend.title = group
  )
  return(p)
}
