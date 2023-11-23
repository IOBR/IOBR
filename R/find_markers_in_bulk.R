


#' Find markers in bulk
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
#' @inheritParams Seurat::RunPCA
#'
#' @return
#' @export
#'
#' @examples
find_markers_in_bulk<-function(pdata, eset, group, id_pdata = "ID", nfeatures = 2000, top_n = 20, thresh.use = 0.25, only.pos = TRUE, min.pct = 0.25, npcs = 50){

  library(Seurat)
  if(dim(pdata)[2] == 2){
    pdata[,"newid"] <- pdata[, id_pdata]
  }

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
  sce <- Seurat:: ScaleData(object = sce,
                   # vars.to.regress = c('nCount_RNA'),
                   model.use = 'linear',
                   use.umi = FALSE)
  sce <- Seurat::FindVariableFeatures(object = sce, nfeatures = nfeatures)
  # length(VariableFeatures(sce))
  if (ncol(sce) < npcs) {
    message("The sample number is less than PC number setting, using the sample number -1 (you can set with npcs)")
    npcs = ncol(sce) - 1
  }
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
