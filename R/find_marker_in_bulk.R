


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
find_markers_in_bulk<-function(pdata, eset, group, id_pdata = "ID", nfeatures = 2000, top_n = 20, thresh.use = 0.25, only.pos = TRUE, min.pct = 0.25, npcs = 30){
  
  library(Seurat)
  if(dim(pdata)[2] == 2){
    pdata[,"newid"] <- pdata[, id_pdata]
  }
  
  rownames(pdata) <- NULL
  pdata<-column_to_rownames(pdata, var = id_pdata)
  pdata<-pdata[rownames(pdata)%in%colnames(eset),]
  # print(mhead(pdata))
  
  feas<-rownames(eset)
  feas<-feature_manipulation(data = eset, feature = feas, is_matrix = TRUE)
  
  eset<-eset[rownames(eset)%in%feas, ]
  eset<- eset[,colnames(eset)%in%rownames(pdata)]
  
  # print(mhead(eset))
  
  sce <- Seurat:: CreateSeuratObject(counts = eset,
                                     meta.data = pdata,
                                     min.cells = 0,
                                     min.features = 0,
                                     project = "sce")
  # print(head(sce@meta.data) )
  sce <- Seurat::NormalizeData(sce)
  sce@assays[["RNA"]]@layers[["data"]]<-sce@assays[["RNA"]]@layers[["counts"]]
  sce <- Seurat:: ScaleData(object = sce,
                           
                            use.umi = FALSE)
  sce <- Seurat::FindVariableFeatures(object = sce, nfeatures = nfeatures)
  # length(VariableFeatures(sce))
  sce <- Seurat::RunPCA(object = sce, pc.genes = VariableFeatures(sce), npcs = npcs)
  # sce <- FindNeighbors(object = sce, dims = 1:20, verbose = FALSE)
  # sce <- FindClusters(object = sce, resolution = 0.5,verbose = FALSE)
  meta <- sce@meta.data
  
  Idents(sce) <- as.character( meta[,group])
  print(table(Seurat::Idents(sce)))
  sce.markers <- Seurat::FindAllMarkers(object = sce, only.pos = only.pos, min.pct = min.pct, thresh.use = thresh.use, logfc.threshold = thresh.use)
  # sce.markers %>% group_by(cluster) %>% top_n(top_n, avg_log2FC)
  topN <- sce.markers %>% group_by(cluster) %>% top_n(top_n, avg_log2FC)
  print(topN)
  res<-list("sce" = sce, "markers" = sce.markers, "top_markers" = topN)
  return(res)
}
