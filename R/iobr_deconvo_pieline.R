

#' Title iobr_deconvo_pieline
#'
#' @param eset 
#' @param project 
#' @param array 
#' @param tumor_type 
#' @param path 
#' @param permutation 
#'
#' @return
#' @export
#'
#' @examples
iobr_deconvo_pieline <- function(eset, project, array, tumor_type, path = "1-TME", permutation = 1000){
  
  
  #######################################
  path <- creat_folder(paste0(path))
  #############################################
  eset<- log2eset(eset_tpm)
  #############################################
  
  # tumor_type<-"stad"
  #########################################
  cibersort<-deconvo_tme(eset = eset,method = "cibersort",arrays = array,perm = 1000 )
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
  
  sig_res<-calculate_sig_score(pdata = NULL,
                               eset = eset,
                               signature = signature_collection,
                               method = "integration",
                               mini_gene_count = 2)
  print(paste0( ">>>>> Signature esitmation was completed: ", project))
  save(sig_res,file = paste0(path$abspath,"2-",project,"-Signature-score-mycollection.RData"))
  ########################################
  sig_go_kegg<-calculate_sig_score(pdata = NULL,
                                   eset = eset,
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