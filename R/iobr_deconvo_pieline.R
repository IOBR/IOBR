

#' Tumor Microenvironment Analysis Pipeline
#'
#' The iobr_deconvo_pipeline function executes a comprehensive TME analysis on a gene expression dataset.
#' This analysis includes TME deconvolution using various computational methods, calculation of signature scores,
#' and integration of these data into a cohesive output. It is designed for in-depth exploration of the
#' microenvironment's role in tumor biology based on gene expression profiles.
#'
#' @param eset A gene expression dataset, typically an expression set object, prepared for TME analysis.
#' @param project A character string specifying the project or analysis name, used for output file naming.
#' @param array A logical indicating whether the data comes from an array platform; influences deconvolution methods.
#' @param tumor_type A character string specifying the tumor type, which can tailor certain analysis aspects.
#' @param path A string indicating the output file path; defaults to "1-TME".
#' @param permutation An integer specifying the number of permutations to use in TME deconvolution methods; default is 1000.
#'
#' @return A comprehensive dataset combining TME cell fractions, signature scores, and integrated TME-signature analysis.
#'         This function saves several files to the specified path, documenting various stages of the analysis.
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data("eset_stad", package = "IOBR")
#' eset <- count2tpm(eset_stad)  # Prepare data
#' # Run the pipeline
#' res <- iobr_deconvo_pipeline(eset = eset, project = "STAD", array = FALSE, tumor_type = "stad", path = "1-TME", permutation = 1000)
iobr_deconvo_pipeline <- function(eset, project, array, tumor_type, path = "1-TME", permutation = 1000){


  #######################################
  path <- creat_folder(paste0(path))
  #############################################
  eset<- log2eset(eset)
  #############################################

  # tumor_type<-"stad"
  #########################################
  cibersort<-deconvo_tme(eset = eset,method = "cibersort",arrays = array, perm = permutation )
  epic     <-deconvo_tme(eset = eset,method = "epic",arrays = array)
  mcp      <-deconvo_tme(eset = eset,method = "mcpcounter")
  xcell    <-deconvo_tme(eset = eset,method = "xcell",arrays = array)
  estimate <-deconvo_tme(eset = eset,method = "estimate")
  timer    <-deconvo_tme(eset = eset,method = "timer",group_list = rep(tumor_type,dim(eset)[2]))
  quantiseq<-deconvo_tme(eset = eset,method = "quantiseq", tumor = TRUE, arrays = array, scale_mrna = TRUE)
  ips      <-deconvo_tme(eset = eset,method = "ips",plot= FALSE)

  tme_combine<-cibersort %>%
    inner_join(.,mcp,by       = "ID") %>%
    inner_join(.,xcell,by     = "ID") %>%
    inner_join(.,epic,by      = "ID") %>%
    inner_join(.,estimate,by  = "ID") %>%
    inner_join(.,quantiseq,by = "ID") %>%
    inner_join(.,timer,by     = "ID") %>%
    inner_join(.,ips,by       = "ID")

  # tme_combine<-tme_combine[,-c(grep(colnames(tme_combine),pattern = "Index"))]
  ############################################
  save(tme_combine,file = paste0(path$abspath,"1-",project,"-TME-Cell-fration.RData"))
  print(paste0( ">>>>> TME cell deconvolution was completed: ", project))
  #######################################

  data(signature_collection, package = "IOBR")
  sig_res<-calculate_sig_score(pdata = NULL,
                               eset = eset,
                               signature = signature_collection,
                               adjust_eset = TRUE,
                               method = "integration",
                               mini_gene_count = 2)
  print(paste0( ">>>>> Signature esitmation was completed: ", project))
  save(sig_res,file = paste0(path$abspath,"2-",project,"-Signature-score-mycollection.RData"))
  ########################################
  sig_go_kegg<-calculate_sig_score(pdata = NULL,
                                   eset = eset,
                                   adjust_eset = TRUE,
                                   signature = c(hallmark,go_bp,go_cc,go_mf,kegg,reactome),
                                   method = "ssgsea",
                                   mini_gene_count = 2)


  save(sig_go_kegg,file = paste0(path$abspath,"3-",project,"-Signature-score-Hallmark-GO-KEGG.RData"))
  print(paste0( ">>>>> HALLMARK GO KEGG REACTOME esitmation was completed: ", project))
  #########################################
  tme_sig_combin<-tme_combine %>%
    inner_join(.,sig_res,by = "ID") %>%
    inner_join(.,sig_go_kegg,by = "ID")
  save(tme_sig_combin,file = paste0(path$abspath,"0-",project,"-Merge-TME-Signature-and-Hallmark-GO-KEGG.RData"))
  #########################################
  return(tme_sig_combin)
}
