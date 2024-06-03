


#' Title tcga_rna_preps
#'
#' @description
#' The tcga_rna_preps function is designed to preprocess TCGA RNA-seq data for downstream analysis. Upon executing the function, it performs various operations based on the specified parameters, including modifying sample types, transforming input data, and annotating genes. The function returns the preprocessed RNA-seq gene expression data in the specified output format.
#'
#' @param eset A matrix or data frame representing RNA-seq gene expression data, typically from The Cancer Genome Atlas (TCGA).
#' @param id_type (optional, default = "ensembl"): Specifies the type of gene identifier to be used, accepting values "ensembl" or "symbol".
#' @param input_type  (optional, default = "log2count"): Indicates the input data type, accepting values "log2count" or "count".
#' @param output  (optional, default = "tumor"): Specifies the output sample type, accepting values "tumor" or "tumor_normal".
#' @param output_type  (optional, default = "tpm"): Indicates the output data type, accepting values "tpm", "log2tpm", or "count".
#' @param annotation (optional, default = TRUE): A logical value determining whether gene annotation should be applied.
#'
#' @return A matrix or data frame of the preprocessed gene expression data.
#' @export
#'
#' @author Dongqiang Zeng
#' @examples
#' # Expression matrix downloaded from UCSC Xena datasets: https://xenabrowser.net/datapages/
#' data("eset_stad", package = "IOBR")
#' eset <- tcga_rna_preps(eset = eset_stad,  id_type = "ensembl", input_type = "count", output = "tumor", output_type = "tpm", annotation = TRUE)
#'
tcga_rna_preps <- function(eset, id_type = c("ensembl", "symbol"), input_type = c("log2count", "count"), output = c("tumor", "tumor_normal"),
                           output_type = c("tpm", "log2tpm", "count"), annotation = TRUE){


  if(id_type=="ensembl")  rownames(eset) <- substring(rownames(eset), 1, 15)

  # Revert back to original format because the data from UCSC was log2(x+1)transformed.
  if(input_type=="log2count")  eset<-(2^eset)-1

  #remove ajacent normal sample
  if(output=="tumor"){

    message(">>== sample type: ")
    message(print(table(substring(colnames(eset), 14,16))))

    message(">>== Ajacent normal sample was removed... ")
    eset<- eset[,!substring(colnames(eset), 14,16)=="11A"]

    message(">>== TCGA barcode reduced to 12 digits... ")
    colnames(eset)<- substring(colnames(eset), 1,12)

    message(">>== Presence of tumour samples with the same barcode? ")
    print(summary(duplicated(colnames(eset))))
    eset<- eset[,!duplicated(colnames(eset))]
  }else{
    message(">>>== Both tumour and paracancerous samples were retained and TCGA barcode was retained to 16 digits")
  }

 if(output_type=="tpm"){

   message(">>== count to TPM... ")
   eset_tpm <- count2tpm(countMat = eset, idType = id_type)
   # print(head(eset))

 }else if(output_type=="log2tpm"){

   message(">>== count to TPM and annotation... ")
   eset_tpm <- count2tpm(countMat = eset, idType = id_type)
   eset_tpm <- log2eset(eset_tpm)
   # print(head(eset))
 }else if(output_type=="count"){
   eset_tpm <- eset
   if(id_type=="ensembl" & annotation ){
     message(">>== Annotation... ")
     data("anno_grch38", package = "IOBR")
     eset_tpm <- anno_eset(eset = eset_tpm, annotation = anno_grch38, probe = "id", symbol = "symbol")
   }
 }

  return(eset_tpm)
}




