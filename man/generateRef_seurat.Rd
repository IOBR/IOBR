% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateRef_seurat.R
\name{generateRef_seurat}
\alias{generateRef_seurat}
\title{generateRef_seurat}
\usage{
generateRef_seurat(
  sce,
  celltype = NULL,
  proportion = NULL,
  assay_deg = "RNA",
  slot_deg = "data",
  adjust_assay = FALSE,
  assay_out = "RNA",
  slot_out = "data",
  verbose = FALSE,
  only.pos = TRUE,
  n_ref_genes = 50,
  logfc.threshold = 0.15,
  test.use = "wilcox"
)
}
\arguments{
\item{sce}{Seurat object containing single-cell RNA-seq data.}

\item{celltype}{(Optional) A character vector specifying cell type for which marker genes will be identified.}

\item{proportion}{(Optional) A numeric value specifying the proportion of cells to be randomly selected for analysis.}

\item{assay_deg}{(Default = "RNA") A character string specifying the assay to be used for finding markers.}

\item{slot_deg}{(Default = "data") A character string specifying the slot to be used for finding markers.}

\item{adjust_assay}{(Default = FALSE) A logical value specifying whether to adjust the assay.}

\item{assay_out}{(Default = "RNA") A character string specifying the assay used in the output.}

\item{slot_out}{(Default = "data") A character string specifying the slot used in the output.}

\item{verbose}{(Default = FALSE) A logical value specifying whether to print verbose messages.}

\item{only.pos}{(Default = TRUE) A logical value specifying whether to consider only positive results.}

\item{n_ref_genes}{(Default = 50) An integer specifying the number of reference genes to be selected.}

\item{logfc.threshold}{(Default = 0.15) A numeric value specifying the threshold for log-fold change.}

\item{test.use}{(Default = "wilcox") A character string specifying the statistical test to be used.}
}
\value{
A matrix containing aggregated expression data for the reference genes.
}
\description{
The "generateRef_seurat" function takes a Seurat object "sce" and additional parameters to perform various operations for generating reference gene expression data. The function allows for specifying cell types, proportions, assays, preprocessing options, and statistical testing parameters.
}
\examples{

pbmc.data <- Read10X(data.dir = "E:/12-pkg-dev/IOBR-Project/8-IOBR2-Vignette/0-data/pbmc/filtered_gene_bc_matrices/hg19")
Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- ScaleData(pbmc, features =  rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc$celltype <- paste0("celltype_", pbmc$seurat_clusters)

# generate reference matrix
sm<- generateRef_seurat(sce = pbmc, celltype = "celltype", slot_out = "data")

#load the bulk-seq data
data(eset_stad, package = "IOBR")
eset <- count2tpm(countMat = eset_stad, source = "local", idType = "ensembl")
svr<-deconvo_tme(eset = eset, reference  = sm,  method = "svr", arrays  = FALSE,absolute.mode = FALSE, perm = 100)

}
\author{
Dongqiang Zeng
}
