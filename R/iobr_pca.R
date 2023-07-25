




#' Title
#'
#' @param data
#' @param is.matrix
#' @param scale
#' @param is.log
#' @param pdata
#' @param id_pdata
#' @param group
#' @param cols
#' @param palette
#' @param repel
#' @param ncp
#' @param axes
#' @param addEllipses
#'
#' @return
#' @export
#'
#' @examples
iobr_pca <- function(data, is.matrix = TRUE, scale = TRUE, is.log = FALSE, pdata, id_pdata = "ID", group = NULL,
                     cols = "normal", palette = "jama", repel = FALSE, ncp = 6, axes = c(1, 2), addEllipses = TRUE){

  if(is.matrix) data <- t(data)
  if(is.log) data <- log2eset(data)
  if(scale) data <- scale(data)

  res.pca <- FactoMineR::PCA(data, ncp = ncp, graph = FALSE)

  colnames(pdata)[which(colnames(pdata)==id_pdata)] <- "id"
  pdata <- pdata[pdata$id%in%rownames(data), ]
  pdata <- pdata[match(rownames(data), pdata$id), ]

  message(print(table(pdata[, group])))
  ##########################################

  cols <- get_cols(cols = cols, palette = palette, show_col = FALSE, seed = 123)

  cols <- cols[1:length(unique(pdata[, group]))]
  #########################################
  p <- factoextra:: fviz_pca_ind(res.pca,
                                axes = axes,
                                geom.ind     = "point", # show points only (nbut not "text")
                                col.ind      = pdata[, group], # color by groups
                                palette      = cols,
                                repel        = repel,
                                addEllipses  = addEllipses, # Concentration ellipses
                                legend.title = group)
  return(p)

}
