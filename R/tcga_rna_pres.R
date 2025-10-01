#' Preprocess TCGA RNA-seq Data
#'
#' Preprocesses TCGA RNA-seq data by modifying sample types, transforming data,
#' and annotating genes based on specified parameters.
#'
#' @param eset RNA-seq gene expression matrix from TCGA.
#' @param id_type Gene identifier type: "ensembl" or "symbol". Default is "ensembl".
#' @param input_type Input data type: "log2count" or "count". Default is "log2count".
#' @param output Sample type: "tumor" or "tumor_normal". Default is "tumor".
#' @param output_type Output data type: "tpm", "log2tpm", or "count". Default is "tpm".
#' @param annotation Logical for gene annotation. Default is TRUE.
#'
#' @return Preprocessed gene expression matrix.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("eset_stad", package = "IOBR")
#' eset <- tcga_rna_pres(eset = eset_stad, id_type = "ensembl", input_type = "count",
#'                       output = "tumor", output_type = "tpm", annotation = TRUE)
#'
tcga_rna_preps <- function(eset, id_type = c("ensembl", "symbol"), input_type = c("log2count", "count"), output = c("tumor", "tumor_normal"),
                           output_type = c("tpm", "log2tpm", "count"), annotation = TRUE) {
  if (id_type == "ensembl") rownames(eset) <- substring(rownames(eset), 1, 15)

  # Revert back to original format because the data from UCSC was log2(x+1)transformed.
  if (input_type == "log2count") eset <- (2^eset) - 1

  # remove ajacent normal sample
  if (output == "tumor") {
    message(">>== sample type: ")
    message(print(table(substring(colnames(eset), 14, 16))))

    message(">>== Ajacent normal sample was removed... ")
    eset <- eset[, !substring(colnames(eset), 14, 16) == "11A"]

    message(">>== TCGA barcode reduced to 12 digits... ")
    colnames(eset) <- substring(colnames(eset), 1, 12)

    message(">>== Presence of tumour samples with the same barcode? ")
    print(summary(duplicated(colnames(eset))))
    eset <- eset[, !duplicated(colnames(eset))]
  } else {
    message(">>>== Both tumour and paracancerous samples were retained and TCGA barcode was retained to 16 digits")
  }

  if (output_type == "tpm") {
    message(">>== count to TPM... ")
    eset_tpm <- count2tpm(countMat = eset, idType = id_type)
    # print(head(eset))
  } else if (output_type == "log2tpm") {
    message(">>== count to TPM and annotation... ")
    eset_tpm <- count2tpm(countMat = eset, idType = id_type)
    eset_tpm <- log2eset(eset_tpm)
    # print(head(eset))
  } else if (output_type == "count") {
    eset_tpm <- eset
    if (id_type == "ensembl" & annotation) {
      message(">>== Annotation... ")
      data("anno_grch38", package = "IOBR")
      eset_tpm <- anno_eset(eset = eset_tpm, annotation = anno_grch38, probe = "id", symbol = "symbol")
    }
  }

  return(eset_tpm)
}
