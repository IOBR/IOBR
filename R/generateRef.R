#' Generate Reference Signature Matrix
#'
#' @description
#' Generates a reference signature matrix for cell types based on
#' differential expression analysis. Supports both limma for normalized data
#' and DESeq2 for raw count data.
#'
#' @param dds Matrix. Raw count data from RNA-seq. Required if
#'   `method = "DESeq2"`.
#' @param pheno Character vector. Cell type class of the samples.
#' @param FDR Numeric. Genes with BH adjusted p-value < FDR are considered
#'   significant. Default is 0.05.
#' @param dat Matrix or data frame. Normalized transcript quantification data
#'   (e.g., FPKM, TPM).
#' @param method Character. Method for differential expression: `"limma"` or
#'   `"DESeq2"`. Default is `"limma"`.
#'
#' @return List containing:
#'   - `reference_matrix`: Data frame of median expression for significant
#'     genes across cell types.
#'   - `G`: Optimal number of probes minimizing condition number.
#'   - `condition_number`: Minimum condition number.
#'   - `whole_matrix`: Full median expression matrix.
#'
#' @export
#'
#' @examples
#' expressionData <- matrix(runif(1000 * 4, min = 0, max = 10), ncol = 4)
#' rownames(expressionData) <- paste("Gene", 1:1000, sep = "_")
#' colnames(expressionData) <- paste("Sample", 1:4, sep = "_")
#'
#' phenotype <- c("celltype1", "celltype2", "celltype1", "celltype2")
#'
#' rawCountData <- matrix(sample(1:100, 1000 * 4, replace = TRUE), ncol = 4)
#' rownames(rawCountData) <- paste("Gene", 1:1000, sep = "_")
#' colnames(rawCountData) <- paste("Sample", 1:4, sep = "_")
#'
#' result <- generateRef(
#'   dds = rawCountData, pheno = phenotype,
#'   FDR = 0.05, dat = expressionData, method = "DESeq2"
#' )
generateRef <- function(dds, pheno, FDR = 0.05, dat, method = "limma") {
  method <- rlang::arg_match(method, c("limma", "DESeq2"))

  cli::cli_alert_info("Running differentially expressed genes using {method}")

  res <- switch(method,
    limma = generateRef_limma(dat, pheno, FDR),
    DESeq2 = generateRef_DEseq2(dat, pheno, FDR, dds)
  )

  ref <- res$reference_matrix
  res$reference_matrix <- ref[, -1, drop = FALSE]

  res
}
