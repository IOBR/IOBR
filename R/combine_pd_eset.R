





#' Title combine_pd_eset
#'
#' @description combine_pd_eset combines the expression set (eset) with phenotype data (pdata).
#' @param eset Expression set matrix or data frame.
#' @param pdata Phenotype data.
#' @param id_pdata String indicating the column name in pdata that corresponds to the eset column.
#' @param feas Vector of feature names to include in the combined dataset. Default is NULL, which includes all features.
#' @param feature_manipulation Logical value indicating whether to perform feature manipulation on the data. Default is TRUE.
#' @param scale Logical value indicating whether to scale the expression set. Default is TRUE.
#' @param choose_who_when_duplicate String indicating the preference when encountering duplicate IDs in pdata and eset.
#'
#' @return
#' @export
#'
#' @examples
#' # Loading TCGA-STAD expresion data(raw count matrix)
#' data("eset_stad", package = "IOBR")
#' # transform count data to tpm
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' colnames(eset) <- substring(colnames(eset), 1, 12)
#' # Loading TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' eset <- count2tpm(eset_stad)
#' colnames(eset) <- substring(colnames(eset), 1, 12)
#' input <- combine_pd_eset(eset = eset, pdata = sig_stad, feas = unique(unlist(signature_collection)))
combine_pd_eset<- function(eset, pdata, id_pdata = "ID", feas = NULL,  feature_manipulation = TRUE, scale = TRUE, choose_who_when_duplicate = "eset"){

  if(!is.null(feas)){
    eset<-eset[rownames(eset)%in%feas, ]
  }

  feas_filter<-feature_manipulation(data = eset, feature = rownames(eset), is_matrix = T)

  eset<-eset[rownames(eset)%in%feas_filter, ]
  eset<-t(eset)
  if(scale) eset<-scale(eset, scale = T, center = T)
  eset<-rownames_to_column(as.data.frame(eset), var = "ID")

  colnames(pdata)[which(colnames(pdata)==id_pdata)]<- "ID"

  if(choose_who_when_duplicate=="eset"){
    choose<-"X"
  }else{
    choose="Y"
  }
  pd_eset<- merge_duplicate(pdata, eset, by.x = "ID", by.y = "ID", choose = choose, all = FALSE)
  return(pd_eset)

}
