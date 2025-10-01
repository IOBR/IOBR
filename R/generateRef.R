#' Generate Reference Signature Matrix
#'
#' This function generates a reference signature matrix for cell types based on differential expression analysis.
#' It can use either 'limma' for normalized data or 'DESeq2' for raw count data. The function identifies genes
#' with significant expression changes at a specified false discovery rate (FDR) threshold.
#'
#' @param dds raw count data from RNA-seq; Necessary if used the method DESeq2
#' @param pheno character vector; cell type class of the samples
#' @param FDR numeric; genes with BH adjust p value < FDR are considered significant.
#' @param dat data frame or matrix; normalized transcript quantification data (like FPKM, TPM). Note: cell's median expression level of the identified probes will be the output of reference_matrix.
#' @param method limma or DESeq2
#'
#' @return A list containing the reference signature matrix and possibly other elements depending on the analysis method used.
#'         The cells of the matrix represent the median expression level of identified significant genes across samples grouped by cell type.
#' @import DESeq2
#' @export
#'
#' @examples
#' # Simulate expression data for 1000 genes across 4 samples
#' expressionData <- matrix(runif(1000 * 4, min = 0, max = 10), ncol = 4)
#' rownames(expressionData) <- paste("Gene", 1:1000, sep = "_")
#' colnames(expressionData) <- paste("Sample", 1:4, sep = "_")
#'
#' # Create phenotype data for the samples
#' phenotype <- c("celltype1", "celltype2", "celltype1", "celltype2")
#'
#' # Simulate raw count data for 1000 genes across 4 samples
#' rawCountData <- matrix(sample(1:100, 1000 * 4, replace = TRUE), ncol = 4)
#' rownames(rawCountData) <- paste("Gene", 1:1000, sep = "_")
#' colnames(rawCountData) <- paste("Sample", 1:4, sep = "_")
#'
#' # Create column data for building a DESeqDataSet
#' library(DESeq2)
#' colData <- data.frame(
#'   celltype = phenotype,
#'   condition = c("treated", "control", "treated", "control"),
#'   row.names = colnames(rawCountData)
#' )
#' # Assuming the design matrix is based on the condition
#' dds_object <- DESeqDataSetFromMatrix(countData = rawCountData, colData = colData, design = ~condition)
generateRef <- function(dds, pheno, FDR = 0.05, dat, method = "limma") {
  print(message(paste0("\n", ">>> Running differentially expressed genes using ", method)))
  res <- switch(method,
    limma = generateRef_limma(dat, pheno, FDR),
    DESeq2 = generateRef_DEseq2(dat, pheno, FDR, dds)
  )
  ref <- res$reference_matrix
  res$reference_matrix <- ref[, -1]
  return(res)
}
