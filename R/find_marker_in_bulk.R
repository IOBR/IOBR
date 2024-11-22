


#' Find markers in bulk
#' 
#' @description The goal of this function is to find relevant results from the given gene expression data and meta information.
#'
#' @param pdata A data frame containing the meta information of the samples.
#' @param eset A matrix containing the gene expression data or signature score.
#' @param group A string representing the column name for grouping.
#' @param id_pdata A string representing the column name for sample IDs, default is "ID".
#' @param nfeatures A numeric value indicating the top n features to select from the variable features, default is 2000.
#' @param top_n A numeric value representing the top n markers to select in each cluster, default is 20.
#' @param thresh.use A numeric value representing the marker selection threshold, default is 0.25.
#' @param only.pos A logical value indicating whether to select only positive markers, default is TRUE.
#' @param min.pct A numeric value representing the minimum percentage threshold for marker selection, default is 0.25.
#'
#' @import Seurat
#' @import dplyr
#' @return A list containing the Seurat object (`sce`), all markers identified (`markers`), and the top markers per group (`top_markers`).
#' @export
#'
#' @examples
#'
#' # loading expression data
#' data("eset_tme_stad", package = "IOBR")
#' colnames(eset_tme_stad) <- substring(colnames(eset_tme_stad), 1, 12)
#'
#' data("pdata_sig_tme", package = "IOBR")
#' res <- find_markers_in_bulk(pdata = pdata_sig_tme, eset = eset_tme_stad, group = "TMEcluster")
#'
#' # extracting top 15 markers of each TME clusters
#' top15 <-  res$top_markers %>% dplyr:: group_by(cluster) %>%  dplyr::top_n(15, avg_log2FC)
#'
#' # visualization
#' cols <- c('#2692a4','#fc0d3a','#ffbe0b')
#' DoHeatmap(res$sce, top15$gene, group.colors = cols )+ scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)))
#'
find_markers_in_bulk <- function(pdata, eset, group, id_pdata = "ID", nfeatures = 2000, top_n = 20, thresh.use = 0.25, only.pos = TRUE, min.pct = 0.25, npcs = 30) {
  library(Seurat)
  library(Matrix)
  
  if (!inherits(eset, "dgCMatrix")) {
    eset <- as(as.matrix(eset), "dgCMatrix")
  }
  
  if (dim(pdata)[2] == 2) {
    pdata[,"newid"] <- pdata[, id_pdata]
  }
  
  rownames(pdata) <- NULL
  pdata <- column_to_rownames(pdata, var = id_pdata)
  pdata <- pdata[rownames(pdata) %in% colnames(eset),]
  
  sce <- Seurat::CreateSeuratObject(counts = eset,
                                    meta.data = pdata,
                                    min.cells = 0,
                                    min.features = 0,
                                    project = "sce")
  
  seurat_version <- packageVersion("Seurat")
  
  if (seurat_version >= "5.0.0") {
    sce <- Seurat::NormalizeData(sce)
    sce@assays[["RNA"]]@layers[["data"]] <- sce@assays[["RNA"]]@layers[["counts"]]
    sce <- Seurat::ScaleData(sce, layer = "data")
  } else {
    sce <- Seurat::NormalizeData(sce)
    sce <- Seurat::ScaleData(sce)
  }
  
  tryCatch({
    sce <- Seurat::FindVariableFeatures(
      object = sce, 
      selection.method = "vst", 
      mean.cutoff = c(0.1, 8), 
      dispersion.cutoff = c(1, Inf)
    )
  }, error = function(e) {
    message("Error in FindVariableFeatures: ", e$message)
    gene_var <- apply(as.matrix(sce@assays[["RNA"]]@counts), 1, var)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:nfeatures]
    VariableFeatures(sce) <- top_genes
  })
  
  sce <- Seurat::RunPCA(object = sce, features = VariableFeatures(sce), npcs = npcs)
  
  meta <- sce@meta.data
  Idents(sce) <- as.character(meta[, group])
  
  sce.markers <- Seurat::FindAllMarkers(object = sce, only.pos = only.pos, min.pct = min.pct, thresh.use = thresh.use, logfc.threshold = thresh.use)
  topN <- sce.markers %>% group_by(cluster) %>% top_n(top_n, avg_log2FC)
  
  res <- list("sce" = sce, "markers" = sce.markers, "top_markers" = topN)
  return(res)
}
