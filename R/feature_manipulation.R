






#' Filtered features with multiple methods
#'
#' @param data data with features in the columns
#' @param feature features with NA or Infinity value
#' @return filtered features
#' @export
#' @author Dongqiang Zeng
#' @examples

feature_manipulation<-function(data = data,feature = feature){
  data<-as.data.frame(data)
  #remove NA variables
  print(paste(">>> Is NA exist: ",  sum(is.na(data[,feature]))))

  if (sum(is.na(data))>0) {
    nn<-as.data.frame(t(data[,feature]))
    delete_vars<- rownames(nn)[!complete.cases(nn)]
    feature<-feature[!feature%in%delete_vars]
  }
  #####################################

  #remove non-numeric variables
  message(paste0(">>>> Is nonnumeric variables exist ? >>>>"))
  print(summary(lapply(data[,feature],function(x) class(x))=="numeric"))
  fea_class<-as.vector(lapply(data[,feature],function(x) class(x))=="numeric")
  feature<-feature[fea_class]

  #remove infinited variables
  message(paste0(">>>> Is -Inf variables exist ? >>>>"))
  print(summary(lapply(data[,feature],function(x) min(x))==-Inf))
  fea_class<-as.vector(lapply(data[,feature],function(x) min(x))==-Inf)
  feature<-feature[!fea_class]

  message(paste0(">>>> Is Inf variables exist ? >>>>"))
  print(summary(lapply(data[,feature],function(x) max(x))== Inf))
  fea_class<-as.vector(lapply(data[,feature],function(x) max(x))== Inf)
  feature<-feature[!fea_class]

  #remove variables with same number
  sd<-apply(data[,feature],2,function(x) sd(x)==0)
  message(paste0(">>> Variables have sd = 0 :  "))
  print(summary(sd))
  feature<-feature[!sd]

  return(feature)
}
