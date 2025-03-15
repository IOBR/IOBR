


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
  library(tibble)
  library(dplyr)

  # 转换输入数据格式
  if (!inherits(eset, "dgCMatrix")) {
    eset <- as(as.matrix(eset), "dgCMatrix")
  }

  # 处理元数据
  if (ncol(pdata) == 2) {
    pdata[,"newid"] <- pdata[, id_pdata]
  }
  pdata <- tibble::column_to_rownames(pdata, var = id_pdata)
  pdata <- pdata[rownames(pdata) %in% colnames(eset), ]

  # 创建Seurat对象
  sce <- CreateSeuratObject(
    counts = eset,
    meta.data = pdata,
    min.cells = 0,
    min.features = 0,
    project = "sce"
  )

  # 获取版本号并进行规范比较
  seurat_version <- packageVersion("Seurat")
  v5_threshold <- package_version("5.0.0")

  # 版本适配处理
  if (seurat_version >= v5_threshold) {
    # Seurat v5+ 的处理逻辑
    message("Using Seurat v5+ workflow")
    sce <- NormalizeData(sce)
    sce[["RNA"]]$data <- LayerData(sce, assay = "RNA", layer = "counts")  # 兼容v5的Layer访问
    sce <- ScaleData(sce, layer = "data")
  } else {
    # Seurat v4及以下版本的处理逻辑
    message("Using Seurat v4 workflow")
    sce <- NormalizeData(sce)
    sce <- ScaleData(sce)
  }

  # 特征选择
  tryCatch({
    sce <- FindVariableFeatures(
      object = sce,
      selection.method = "vst",
      nfeatures = nfeatures,
      mean.cutoff = c(0.1, 8),
      dispersion.cutoff = c(1, Inf)
    )
  }, error = function(e) {
    message("Error in FindVariableFeatures: ", e$message)
    counts <- GetAssayData(sce, assay = "RNA", slot = "counts")
    gene_var <- Matrix::rowVars(counts)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:nfeatures]
    VariableFeatures(sce) <<- top_genes  # 使用 <<- 修改父环境变量
  })

  # 降维和标记物识别
  sce <- RunPCA(object = sce, features = VariableFeatures(sce), npcs = npcs)
  Idents(sce) <- as.character(sce[[group, drop = TRUE]])

  # 查找差异表达基因
  sce.markers <- FindAllMarkers(
    object = sce,
    only.pos = only.pos,
    min.pct = min.pct,
    logfc.threshold = thresh.use,
    test.use = "wilcox"
  )

  # 提取top标记物
  topN <- sce.markers %>%
    group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = top_n)

  return(list(sce = sce, markers = sce.markers, top_markers = topN))
}

