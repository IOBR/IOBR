


#' Title
#'
#' @param pdata
#' @param eset
#' @param group
#' @param id_pdata
#' @param nfeatures
#' @param top_n
#' @param thresh.use
#' @param only.pos
#' @param min.pct
#'
#' @return
#' @export
#'
#' @examples
find_markers_in_bulk<-function(pdata, eset, group, id_pdata = "ID", nfeatures = 200, top_n = 20, thresh.use = 0.25, only.pos = TRUE, min.pct = 0.25){


  pdata<-column_to_rownames(pdata, var = id_pdata)
  pdata<-pdata[rownames(pdata)%in%colnames(eset),]

  feas<-rownames(eset)
  feas<-feature_manipulation(data = eset, feature = feas, is_matrix = T)

  eset<-eset[rownames(eset)%in%feas, ]
  eset<- eset[,colnames(eset)%in%rownames(pdata)]
  sce <- CreateSeuratObject(counts = eset,
                            meta.data = pdata,
                            min.cells = 0,
                            min.features = 0,
                            project = "sce")
  # print(head(sce@meta.data) )
  sce <- ScaleData(object = sce,
                   # vars.to.regress = c('nCount_RNA'),
                   model.use = 'linear',
                   use.umi = FALSE)
  sce <- FindVariableFeatures(object = sce, nfeatures = nfeatures)
  # length(VariableFeatures(sce))
  sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
  # sce <- FindNeighbors(object = sce, dims = 1:20, verbose = FALSE)
  # sce <- FindClusters(object = sce, resolution = 0.5,verbose = FALSE)
  meta<-sce@meta.data
  Idents(sce)<-as.character( meta[,group])
  print(table(Idents(sce)))
  sce.markers <- FindAllMarkers(object = sce, only.pos = only.pos, min.pct = min.pct, thresh.use = thresh.use, logfc.threshold = thresh.use)
  # sce.markers %>% group_by(cluster) %>% top_n(top_n, avg_log2FC)
  topN <- sce.markers %>% group_by(cluster) %>% top_n(top_n, avg_log2FC)
  print(topN)
  res<-list("sce" = sce, "markers" = sce.markers, "top_markers" = topN)
  return(res)
}
