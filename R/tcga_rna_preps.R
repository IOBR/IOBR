#' Preprocess TCGA RNA-seq Data
#'
#' @description
#' Preprocesses TCGA RNA-seq data by modifying sample types, transforming data,
#' and annotating genes based on specified parameters.
#'
#' @param eset Matrix or data frame. RNA-seq gene expression matrix from TCGA.
#' @param id_type Character. Gene identifier type: "ensembl" or "symbol".
#'   Default is "ensembl".
#' @param input_type Character. Input data type: "log2count" or "count".
#'   Default is "log2count".
#' @param output Character. Sample type: "tumor" or "tumor_normal".
#'   Default is "tumor".
#' @param output_type Character. Output data type: "tpm", "log2tpm", or "count".
#'   Default is "tpm".
#' @param annotation Logical. Whether to perform gene annotation. Default is TRUE.
#'
#' @return Preprocessed gene expression matrix.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' eset_stad <- load_data("eset_stad")
#' eset <- tcga_rna_preps(
#'   eset = eset_stad, id_type = "ensembl", input_type = "count",
#'   output = "tumor", output_type = "tpm", annotation = TRUE
#' )
tcga_rna_preps <- function(eset,
                           id_type = c("ensembl", "symbol"),
                           input_type = c("log2count", "count"),
                           output = c("tumor", "tumor_normal"),
                           output_type = c("tpm", "log2tpm", "count"),
                           annotation = TRUE) {
  id_type <- rlang::arg_match(id_type)
  input_type <- rlang::arg_match(input_type)
  output <- rlang::arg_match(output)
  output_type <- rlang::arg_match(output_type)

  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort("{.arg eset} must be a matrix or data frame.")
  }

  if (id_type == "ensembl") {
    rownames(eset) <- substring(rownames(eset), 1, 15)
  }

  if (input_type == "log2count") {
    eset <- (2^eset) - 1
  }

  if (output == "tumor") {
    cli::cli_alert_info("Sample type distribution:")
    message(paste(capture.output(table(substring(colnames(eset), 14, 16))), collapse = "\n"))

    cli::cli_alert_info("Adjacent normal samples removed")
    eset <- eset[, !substring(colnames(eset), 14, 16) == "11A"]

    cli::cli_alert_info("TCGA barcode reduced to 12 digits")
    colnames(eset) <- substring(colnames(eset), 1, 12)

    cli::cli_alert_info("Duplicate barcode check:")
    message(paste(capture.output(summary(duplicated(colnames(eset)))), collapse = "\n"))
    eset <- eset[, !duplicated(colnames(eset))]
  } else {
    cli::cli_alert_info(
      "Both tumor and paracancerous samples retained; TCGA barcode kept to 16 digits"
    )
  }

  if (output_type == "tpm") {
    cli::cli_alert_info("Converting count to TPM")
    eset_tpm <- count2tpm(countMat = eset, idType = id_type)
  } else if (output_type == "log2tpm") {
    cli::cli_alert_info("Converting count to TPM and applying log2 transformation")
    eset_tpm <- count2tpm(countMat = eset, idType = id_type)
    eset_tpm <- log2eset(eset_tpm)
  } else {
    eset_tpm <- eset
    if (id_type == "ensembl" && annotation) {
      cli::cli_alert_info("Annotating genes")
      anno_grch38 <- load_data("anno_grch38")
      eset_tpm <- anno_eset(
        eset = eset_tpm,
        annotation = anno_grch38,
        probe = "id",
        symbol = "symbol"
      )
    }
  }

  eset_tpm
}
