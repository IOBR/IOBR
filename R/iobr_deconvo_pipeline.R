#' Tumor Microenvironment (TME) Deconvolution Pipeline
#'
#' Executes an integrated TME analysis on a gene expression matrix: performs
#' immune/stromal cell deconvolution using multiple algorithms, computes
#' signature scores, and aggregates results. Designed for exploratory
#' immunogenomic profiling.
#'
#' @param eset Numeric matrix. Gene expression (TPM/log scale) with genes in rows.
#' @param project Character. Project name (used in output naming).
#' @param array Logical. Whether data originated from an array platform. Affects deconvolution choices.
#' @param tumor_type Character. Tumor type code (e.g., "stad") used by certain methods.
#' @param path Character. Output directory. Default is NULL (uses tempdir()).
#' @param permutation Integer. Number of permutations for CIBERSORT (and similar). Default is 1000.
#'
#' @return Data frame integrating cell fractions and signature scores (also writes intermediate outputs to disk).
#' @author Dongqiang Zeng
#' @export
#' @examples
#' lm22 <- load_data("lm22")
#' cancer_genes <- load_data("cancer_type_genes")
#' if (!is.null(lm22) && !is.null(cancer_genes)) {
#'   set.seed(123)
#'   # Create simulated data using lm22 genes to avoid deconvolution failures
#'   genes <- rownames(lm22)
#'   # Add xcell genes so xCell doesn't crash from insufficient genes
#'   xcell <- load_data("xCell.data")
#'   if (!is.null(xcell)) genes <- unique(c(genes, xcell$genes))
#'   # Add cancer type genes for TIMER
#'   genes <- unique(c(genes, cancer_genes[["stad"]]))
#'   
#'   eset <- matrix(runif(length(genes) * 2), nrow = length(genes), ncol = 2)
#'   rownames(eset) <- genes
#'   colnames(eset) <- paste0("Sample", 1:2)
#'   
#'   res <- iobr_deconvo_pipeline(
#'     eset = eset, project = "TEST",
#'     array = FALSE, tumor_type = "stad",
#'     path = tempdir(), permutation = 2
#'   )
#'   if (!is.null(res)) head(res)
#' }
#'
iobr_deconvo_pipeline <- function(eset, project, array, tumor_type, path = NULL, permutation = 1000) {

  if (is.null(eset)) return(NULL)

  if (is.null(path)) {
    path <- tempdir()
  }
  #######################################
  path <- creat_folder(paste0(path))
  #######################################
  eset <- log2eset(eset)
  #######################################

  # tumor_type<-"stad"
  #########################################
  cibersort <- deconvo_tme(eset = eset, method = "cibersort", arrays = array, perm = permutation)
  epic <- deconvo_tme(eset = eset, method = "epic", arrays = array)
  mcp <- deconvo_tme(eset = eset, method = "mcpcounter")
  xcell <- deconvo_tme(eset = eset, method = "xcell", arrays = array)
  estimate <- deconvo_tme(eset = eset, method = "estimate")
  timer <- deconvo_tme(eset = eset, method = "timer", group_list = rep(tumor_type, dim(eset)[2]))
  quantiseq <- deconvo_tme(eset = eset, method = "quantiseq", tumor = TRUE, arrays = array, scale_mrna = TRUE)
  ips <- deconvo_tme(eset = eset, method = "ips", plot = FALSE)

  res_list <- list(cibersort, mcp, xcell, epic, estimate, quantiseq, timer, ips)
  res_list <- res_list[!sapply(res_list, is.null)]

  if (length(res_list) == 0) {
    cli::cli_alert_warning("All deconvolution methods returned NULL (likely due to no internet connection).")
    return(NULL)
  }

  tme_combine <- res_list[[1]]
  if (length(res_list) > 1) {
    for (i in 2:length(res_list)) {
      tme_combine <- dplyr::inner_join(tme_combine, res_list[[i]], by = "ID")
    }
  }

  # tme_combine<-tme_combine[,-c(grep(colnames(tme_combine),pattern = "Index"))]
  ############################################
  save(tme_combine, file = paste0(path$abspath, "1-", project, "-TME-Cell-fration.RData"))
  print(paste0(">>>>> TME cell deconvolution was completed: ", project))
  #######################################

  signature_collection <- load_data("signature_collection")
  if (is.null(signature_collection)) return(NULL)

  sig_res <- calculate_sig_score(
    pdata = NULL,
    eset = eset,
    signature = signature_collection,
    adjust_eset = TRUE,
    method = "integration",
    mini_gene_count = 2
  )
  print(paste0(">>>>> Signature esitmation was completed: ", project))
  save(sig_res, file = paste0(path$abspath, "2-", project, "-Signature-score-mycollection.RData"))
  ########################################
  hallmark_data <- load_data("hallmark")
  if (is.null(hallmark_data)) return(NULL)

  go_bp_data <- load_data("go_bp")
  if (is.null(go_bp_data)) return(NULL)

  go_cc_data <- load_data("go_cc")
  if (is.null(go_cc_data)) return(NULL)

  go_mf_data <- load_data("go_mf")
  if (is.null(go_mf_data)) return(NULL)

  kegg_data <- load_data("kegg")
  if (is.null(kegg_data)) return(NULL)

  reactome_data <- load_data("reactome")
  if (is.null(reactome_data)) return(NULL)

  sig_go_kegg <- calculate_sig_score(
    pdata = NULL,
    eset = eset,
    adjust_eset = TRUE,
    signature = c(hallmark_data, go_bp_data, go_cc_data, go_mf_data, kegg_data, reactome_data),
    method = "ssgsea",
    mini_gene_count = 2
  )


  save(sig_go_kegg, file = paste0(path$abspath, "3-", project, "-Signature-score-Hallmark-GO-KEGG.RData"))
  print(paste0(">>>>> HALLMARK GO KEGG REACTOME esitmation was completed: ", project))
  #########################################
  tme_sig_combin <- tme_combine
  if (!is.null(sig_res)) {
    tme_sig_combin <- dplyr::inner_join(tme_sig_combin, sig_res, by = "ID")
  }
  if (!is.null(sig_go_kegg)) {
    tme_sig_combin <- dplyr::inner_join(tme_sig_combin, sig_go_kegg, by = "ID")
  }
  save(tme_sig_combin, file = paste0(path$abspath, "0-", project, "-Merge-TME-Signature-and-Hallmark-GO-KEGG.RData"))
  #########################################
  return(tme_sig_combin)
}
