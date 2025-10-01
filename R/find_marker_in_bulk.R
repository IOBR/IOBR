#' Identify Marker Features in Bulk Expression Data
#'
#' Identifies informative marker features across groups from bulk gene expression or signature score matrices using Seurat workflows. Performs feature selection, scaling, PCA, clustering, and marker discovery.
#'
#' @param pdata Data frame. Sample metadata.
#' @param eset Matrix. Gene expression or signature score matrix.
#' @param group Character. Column name in pdata specifying grouping variable.
#' @param id_pdata Character. Column name for sample IDs. Default is "ID".
#' @param nfeatures Integer. Number of top variable features to select. Default is 2000.
#' @param top_n Integer. Number of top markers to retain per cluster. Default is 20.
#' @param thresh.use Numeric. Threshold for marker selection. Default is 0.25.
#' @param only.pos Logical. Whether to retain only positive markers. Default is TRUE.
#' @param min.pct Numeric. Minimum expression percentage threshold. Default is 0.25.
#' @param npcs Integer. Number of principal components to use. Default is 30.
#'
#' @import Seurat
#' @import dplyr
#' @return List with components: `sce` (Seurat object), `markers` (all markers), `top_markers` (top markers per group).
#' @export
#'
#' @examples
#' data("eset_tme_stad", package = "IOBR")
#' colnames(eset_tme_stad) <- substring(colnames(eset_tme_stad), 1, 12)
#' data("pdata_sig_tme", package = "IOBR")
#' res <- find_markers_in_bulk(pdata = pdata_sig_tme, eset = eset_tme_stad, group = "TMEcluster")
#' # Extract top 15 markers per cluster
#' top15 <- res$top_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(15, avg_log2FC)
#'
find_markers_in_bulk <- function(pdata, eset, group, id_pdata = "ID", nfeatures = 2000, top_n = 20, thresh.use = 0.25, only.pos = TRUE, min.pct = 0.25, npcs = 30) {
  # Check required packages
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed.")
  }

  # 转换输入数据格式
  if (!inherits(eset, "dgCMatrix")) {
    eset <- as(as.matrix(eset), "dgCMatrix")
  }

  # 处理元数据
  if (ncol(pdata) == 2) {
    pdata[, "newid"] <- pdata[, id_pdata]
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
    sce[["RNA"]]$data <- LayerData(sce, assay = "RNA", layer = "counts") # 兼容v5的Layer访问
    sce <- ScaleData(sce, layer = "data")
  } else {
    # Seurat v4及以下版本的处理逻辑
    message("Using Seurat v4 workflow")
    sce <- NormalizeData(sce)
    sce <- ScaleData(sce)
  }

  # 特征选择
  tryCatch(
    {
      sce <- FindVariableFeatures(
        object = sce,
        selection.method = "vst",
        nfeatures = nfeatures,
        mean.cutoff = c(0.1, 8),
        dispersion.cutoff = c(1, Inf)
      )
    },
    error = function(e) {
      message("Error in FindVariableFeatures: ", e$message)
      counts <- GetAssayData(sce, assay = "RNA", slot = "counts")
      gene_var <- Matrix::rowVars(counts)
      top_genes <- names(sort(gene_var, decreasing = TRUE))[1:nfeatures]
      VariableFeatures(sce) <- top_genes # Assign variable features
    }
  )

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
