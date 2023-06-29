





#' Title
#'
#' @param eset expression set
#' @param pdata phenotype data
#' @param id_pdata identifier of pdata which could be matched with column name of eset
#' @param feas selected features of eset to combine with pdata
#' @param feature_manipulation default is T, data with outliers will be filtered
#' @param scale default is true, scaling data
#' @param choose_who_when_duplicate when duplicated variables exist, preserve data of pdata or eset.
#'
#' @return
#' @export
#'
#' @examples
combine_pd_eset<- function(eset, pdata, id_pdata = "ID", feas = NULL,  feature_manipulation = T, scale = T, choose_who_when_duplicate = "eset"){

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
