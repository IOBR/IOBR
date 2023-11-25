#' Generate Reference signature matrix
#'
#' @param dds raw count data from RNA-seq; Necessary if used the method DESeq2
#' @param pheno character vector; cell type class of the samples
#' @param FDR numeric; genes with BH adjust p value < FDR are considered significant.
#' @param dat data frame or matrix; normalized transcript quantification data (like FPKM, TPM). Note: cell's median expression level of the identified probes will be the output of reference_matrix.
#' @param method limma or DESeq2
#'
#' @return
#' @export
#'
#' @examples
generateRef <- function(dds, pheno, FDR = 0.05, dat, method = "limma"){
  print(message(paste0("\n", ">>> Running differentially expressed genes using ", method)))
  res = switch(method,
               limma = generateRef_limma(dat, pheno, FDR),
               DESeq2 = generateRef_DEseq2(dat, pheno, FDR, dds))
  ref<-res$reference_matrix
  res$reference_matrix<-ref[,-1]
  return(res)
}

