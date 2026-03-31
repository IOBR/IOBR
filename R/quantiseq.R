#' Use quanTIseq to Deconvolute a Gene Expression Matrix
#'
#' @description
#' Deconvolutes gene expression data to estimate immune cell fractions using
#' the quanTIseq method. Source code from https://github.com/FFinotello/quanTIseq.
#'
#' @references
#' F. Finotello, C. Mayer, C. Plattner, G. Laschober, D. Rieder,
#' H. Hackl, A. Krogsdam, W. Posch, D. Wilflingseder, S. Sopper, M. Jsselsteijn,
#' D. Johnsons, Y. Xu, Y. Wang, M. E. Sanders, M. V. Estrada, P.
#' Ericsson-Gonzalez, J. Balko, N. F. de Miranda, Z. Trajanoski.
#' "quanTIseq: quantifying immune contexture of human tumors".
#' bioRxiv 223180. https://doi.org/10.1101/223180.
#'
#' @param mix.mat Data frame or matrix. Gene expression matrix with gene symbols
#'   on the first column and sample IDs on the first row. Expression data must
#'   be on non-log scale (TPM for RNA-seq or expression values for microarrays).
#' @param arrays Logical. Whether expression data are from microarrays.
#'   Default is FALSE. If TRUE, the rmgenes parameter is set to "none".
#' @param signame Character. Name of the signature matrix. Currently only
#'   "TIL10" is available. Default is "TIL10".
#' @param tumor Logical. Whether expression data are from tumor samples.
#'   If TRUE, signature genes with high expression in tumor samples are removed.
#'   Default is FALSE.
#' @param mRNAscale Logical. Whether cell fractions must be scaled to account
#'   for cell-type-specific mRNA content. Default is TRUE.
#' @param method Character. Deconvolution method: "hampel", "huber", "bisquare"
#'   for robust regression, or "lsei" for constrained least squares.
#'   Default is "lsei".
#' @param rmgenes Character. Genes to remove: "unassigned" (default), "default",
#'   "none", or "path".
#'
#' @return Data frame with cell fractions for each sample.
#'
#' @import preprocessCore
#' @export
#' @author Finotello F, et al. (adapted for IOBR)
#'
#' @examples
#' \dontrun{
#' # Example usage with TPM data
#' results <- deconvolute_quantiseq.default(mix.mat = tpm_matrix)
#' }
deconvolute_quantiseq.default <- function(mix.mat,
                                          arrays = FALSE,
                                          signame = "TIL10",
                                          tumor = FALSE,
                                          mRNAscale = TRUE,
                                          method = c("lsei", "hampel", "huber", "bisquare"),
                                          rmgenes = "unassigned") {
  
  method <- rlang::arg_match(method)
  
  if (is.null(mix.mat) || nrow(mix.mat) == 0) {
    cli::cli_abort("{.arg mix.mat} cannot be NULL or empty.")
  }
  
  cli::cli_alert_info("Running quanTIseq deconvolution module")
  
  if (rmgenes == "unassigned" && arrays) {
    rmgenes <- "none"
  } else if (rmgenes == "unassigned" && !arrays) {
    rmgenes <- "default"
  }
  
  listsig <- c("TIL10")
  
  if (signame %in% listsig) {
    sig.mat <- quantiseq_data$TIL10_signature
    mRNA <- quantiseq_data$TIL10_mRNA_scaling
    lrmgenes <- quantiseq_data$TIL10_rmgenes
  } else {
    sig.mat.file <- paste0(signame, "_signature.txt")
    sig.mat <- read.table(sig.mat.file, header = TRUE, sep = "\t", row.names = 1)
    
    mRNA.file <- paste0(signame, "_mRNA_scaling.txt")
    mRNA <- read.table(mRNA.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  }
  
  if (!is.numeric(mix.mat[[1, 1]])) {
    cli::cli_abort(c(
      "Wrong input format for the mixture matrix!",
      "i" = "Please follow the instructions in the documentation."
    ))
  }
  
  if (mRNAscale) {
    colnames(mRNA) <- c("celltype", "scaling")
    mRNA <- as.vector(as.matrix(mRNA$scaling[match(colnames(sig.mat), mRNA$celltype)]))
  } else {
    mRNA <- rep(1, ncol(sig.mat))
  }
  
  cli::cli_alert_info("Gene expression normalization and re-annotation (arrays: {arrays})")
  mix.mat <- fixMixture(mix.mat, arrays = arrays)
  
  if (rmgenes != "none") {
    if (signame %in% listsig) {
      n1 <- nrow(sig.mat)
      sig.mat <- sig.mat[!rownames(sig.mat) %in% lrmgenes, , drop = FALSE]
      n2 <- nrow(sig.mat)
      cli::cli_alert_info("Removing {n1 - n2} noisy genes")
    }
  }
  
  if (tumor) {
    if (signame %in% listsig) {
      abgenes <- quantiseq_data$TIL10_TCGA_aberrant_immune_genes
      n1 <- nrow(sig.mat)
      sig.mat <- sig.mat[!rownames(sig.mat) %in% abgenes, , drop = FALSE]
      n2 <- nrow(sig.mat)
      cli::cli_alert_info("Removing {n1 - n2} genes with high expression in tumors")
    }
  }
  
  ns <- nrow(sig.mat)
  us <- length(intersect(rownames(sig.mat), rownames(mix.mat)))
  perc <- round(us * 100 / ns, digits = 2)
  cli::cli_alert_info("Signature genes found in data set: {us}/{ns} ({perc}%)")
  
  cli::cli_alert_info("Mixture deconvolution (method: {method})")
  results1 <- quanTIseq(sig.mat, mix.mat, scaling = mRNA, method = method)
  
  if ("Tregs" %in% colnames(sig.mat) && "T.cells.CD4" %in% colnames(sig.mat) &&
      method == "lsei") {
    minTregs <- 0.02
    i <- which(colnames(sig.mat) == "T.cells.CD4")
    results2 <- quanTIseq(sig.mat[, -i], mix.mat, scaling = mRNA[-i], method = method)
    
    ind <- which(results1[, "Tregs"] < minTregs)
    if (length(ind) > 0) {
      results1[ind, "Tregs"] <- (results2[ind, "Tregs"] + results1[ind, "Tregs"]) / 2
      results1[ind, "T.cells.CD4"] <- pmax(
        0, results1[ind, "T.cells.CD4"] - (results2[ind, "Tregs"] + results1[ind, "Tregs"]) / 2
      )
    }
  }
  
  results <- results1
  results <- results / apply(results, 1, sum)
  
  cli::cli_alert_success("Deconvolution successful!")
  
  results <- data.frame(results)
  results <- cbind(rownames(results), results)
  colnames(results)[1] <- "Sample"
  
  results
}
