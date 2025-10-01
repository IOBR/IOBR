#' Tumor Microenvironment (TME) Deconvolution Pipeline
#'
#' Executes an integrated TME analysis on a gene expression matrix: performs immune/stromal cell deconvolution using multiple algorithms, computes signature scores, and aggregates results. Designed for exploratory immunogenomic profiling.
#'
#' @param eset Numeric matrix. Gene expression (TPM/log scale) with genes in rows.
#' @param project Character. Project name (used in output naming).
#' @param array Logical. Whether data originated from an array platform. Affects deconvolution choices.
#' @param tumor_type Character. Tumor type code (e.g., "stad") used by certain methods.
#' @param path Character. Output directory. Default is "1-TME".
#' @param permutation Integer. Number of permutations for CIBERSORT (and similar). Default is 1000.
#'
#' @return Data frame integrating cell fractions and signature scores (also writes intermediate outputs to disk).
#' @author Dongqiang Zeng
#' @export
#' @examples
#' data("eset_stad", package = "IOBR")
#' eset <- count2tpm(eset_stad)
#' res <- iobr_deconvo_pipeline(eset = eset, project = "STAD", array = FALSE, tumor_type = "stad", path = "1-TME", permutation = 1000)
iobr_deconvo_pipeline <- function(eset, project, array, tumor_type, path = "1-TME", permutation = 1000) {
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

  tme_combine <- cibersort %>%
    inner_join(., mcp, by = "ID") %>%
    inner_join(., xcell, by = "ID") %>%
    inner_join(., epic, by = "ID") %>%
    inner_join(., estimate, by = "ID") %>%
    inner_join(., quantiseq, by = "ID") %>%
    inner_join(., timer, by = "ID") %>%
    inner_join(., ips, by = "ID")

  # tme_combine<-tme_combine[,-c(grep(colnames(tme_combine),pattern = "Index"))]
  ############################################
  save(tme_combine, file = paste0(path$abspath, "1-", project, "-TME-Cell-fration.RData"))
  print(paste0(">>>>> TME cell deconvolution was completed: ", project))
  #######################################

  data(signature_collection, package = "IOBR")
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
  sig_go_kegg <- calculate_sig_score(
    pdata = NULL,
    eset = eset,
    adjust_eset = TRUE,
    signature = c(hallmark, go_bp, go_cc, go_mf, kegg, reactome),
    method = "ssgsea",
    mini_gene_count = 2
  )


  save(sig_go_kegg, file = paste0(path$abspath, "3-", project, "-Signature-score-Hallmark-GO-KEGG.RData"))
  print(paste0(">>>>> HALLMARK GO KEGG REACTOME esitmation was completed: ", project))
  #########################################
  tme_sig_combin <- tme_combine %>%
    inner_join(., sig_res, by = "ID") %>%
    inner_join(., sig_go_kegg, by = "ID")
  save(tme_sig_combin, file = paste0(path$abspath, "0-", project, "-Merge-TME-Signature-and-Hallmark-GO-KEGG.RData"))
  #########################################
  return(tme_sig_combin)
}
