




#' generateRef_seurat
#'
#' @description
#' The "generateRef_seurat" function takes a Seurat object "sce" and additional parameters to perform various operations for generating reference gene expression data. The function allows for specifying cell types, proportions, assays, preprocessing options, and statistical testing parameters.
#'
#' @param sce  Seurat object containing single-cell RNA-seq data.
#' @param proportion  (Optional) A numeric value specifying the proportion of cells to be randomly selected for analysis.
#' @param celltype (Optional) A character vector specifying cell type for which marker genes will be identified.
#' @param adjust_assay (Default = FALSE) A logical value specifying whether to adjust the assay.
#' @param verbose  (Default = FALSE) A logical value specifying whether to print verbose messages.
#' @param only.pos (Default = TRUE) A logical value specifying whether to consider only positive results.
#' @param n_ref_genes  (Default = 50) An integer specifying the number of reference genes to be selected.
#' @param logfc.threshold (Default = 0.15) A numeric value specifying the threshold for log-fold change.
#' @param test.use  (Default = "wilcox") A character string specifying the statistical test to be used.
#' @param assay_deg (Default = "RNA") A character string specifying the assay to be used for finding markers.
#' @param slot_deg (Default = "data") A character string specifying the slot to be used for finding markers.
#' @param assay_out  (Default = "RNA") A character string specifying the assay used in the output.
#' @param slot_out  (Default = "data") A character string specifying the slot used in the output.
#'
#' @author Dongqiang Zeng
#' @return A matrix containing aggregated expression data for the reference genes.
#' @export
#'
#' @examples
#' \donttest{
#' # Load the PBMC dataset
#' pbmc.data <- Read10X(data.dir = "E:/12-pkg-dev/IOBR-Project/8-IOBR2-Vignette/0-data/pbmc/filtered_gene_bc_matrices/hg19")
#' # Initialize the Seurat object with the raw (non-normalized data).
#' pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
#' pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#' pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#' pbmc <- ScaleData(pbmc, features =  rownames(pbmc))
#' pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#' pbmc <- FindNeighbors(pbmc, dims = 1:10)
#' pbmc <- FindClusters(pbmc, resolution = 0.5)
#' pbmc$celltype <- paste0("celltype_", pbmc$seurat_clusters)
#'
#' # generate reference matrix
#' sm<- generateRef_seurat(sce = pbmc, celltype = "celltype", slot_out = "data")
#'
#' #load the bulk-seq data
#' data(eset_stad, package = "IOBR")
#' eset <- count2tpm(countMat = eset_stad, source = "local", idType = "ensembl")
#' svr<-deconvo_tme(eset = eset, reference  = sm,  method = "svr", arrays  = FALSE,absolute.mode = FALSE, perm = 100)
#' }
#'

generateRef_seurat <- function(sce, celltype = NULL, proportion = NULL, assay_deg = "RNA", slot_deg = "data", adjust_assay = FALSE,
                               assay_out = "RNA", slot_out = "data", verbose = FALSE, only.pos = TRUE, n_ref_genes = 50,
                               logfc.threshold = 0.15, test.use = "wilcox"){

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  # if(!is.null(path)){
  #   file_store<-path
  # }else{
  #   file_store<-paste0("generateRef-result")
  # }
  #
  # path<-creat_folder(file_store)

  # if(!file.exists(file_store)) dir.create(file_store)
  # abspath<-paste(getwd(),"/",file_store,"/",sep ="" )

  cat(crayon::green(">>>---Assay used to find markers: \n"))

  if(!is.null(assay_deg)){
    print(paste0(">>>>> ",assay_deg))
  }else{
    print(paste0(">>>>> ", DefaultAssay(sce)))
    DefaultAssay(sce)<- assay_deg
  }


  # find marker genes of clusters------------------------------------------------------
  if(!is.null(celltype)){
    message(paste0(">>> Idents of Seurat object is: ", celltype))
    Idents(sce) <- celltype
    print(table(as.character(Idents(sce))))
    group2<-celltype
  }else{
    group2<-"default"

    cat(crayon::green(">>> Idents of Seurat object is: Default... \n"))
    # message(paste0(">>> Idents of Seurat object is: Default... \n"))
    print(table(as.character(Idents(sce))))
  }


  ##################################
  # help("PrepSCTFindMarkers")

  # if(tolower(assay)=="sct"&&adjust_assay) sce<- PrepSCTFindMarkers(sce)


  ################################################
  if(!is.null(proportion)){


    input<-random_strata_cells(input= sce@meta.data, group = celltype, propotion = proportion)
    cell_input<-rownames(input)

    message(paste0(">>> Subseting cells in each cluster randomly: ",proportion*100, "%... "))
    message(paste0(">>> Final training data: "))
    print(table(as.factor(input[, celltype])))

    # print(cell_input[1:10])
    # subset cells
    sce<-subset(sce, cells = as.character(cell_input))

  }
  ################################################
  #remove features with large name
  ##############################################
  # help("PrepSCTFindMarkers")

  if(tolower(assay_deg)=="sct"&&adjust_assay) sce<- PrepSCTFindMarkers(sce)


  ################################################
  cat(crayon::green(">>> Find markers of each celltype... \n"))
  # help("FindAllMarkers")
  ###################################
  sce.markers <- FindAllMarkers(object          = sce,
                                slot            = slot_deg,
                                assay           = assay_deg,
                                features        = NULL,
                                only.pos        = only.pos,
                                min.pct         = 0.25,
                                thresh.use      = 0.25,
                                verbose         = verbose,
                                logfc.threshold = logfc.threshold,
                                test.use        = test.use,
                                fc.name         = "avg_log2FC")

  # DT::datatable(sce.markers)
  print(head(sce.markers))
  ####################################
  refgene <- sce.markers %>% dplyr:: group_by(cluster) %>%  dplyr::top_n(n_ref_genes, avg_log2FC)

  print(refgene)
  ####################################


  cat(crayon::green(">>>-- Aggreating scRNAseq data...\n"))
  if(is.null(slot_out)) slot<-"counts"

  cat(crayon::green(">>>-- `orig.ident` was set as group. User can define through parameter `celltype` ...\n"))
  bulk <- Seurat:::PseudobulkExpression(object    = sce,
                                        pb.method = 'aggregate',
                                        slot      = slot_out,
                                        assays    = assay_out,
                                        group.by  = celltype)
  bulk<-bulk[[1]]

  # print(bulk)
  bulk <- bulk[rownames(bulk)%in%refgene$gene, ]

  return(bulk)

}


