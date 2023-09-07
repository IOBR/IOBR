
#' feature_manipulation - Manipulate features in a matrix or dataframe
#'
#' @description This function allows for manipulation of features in a given matrix or dataframe. It can remove variables with missing values, non-numeric variables, infinite variables, and variables with zero standard deviation.
#' @param data The matrix or dataframe containing the features to be manipulated.
#' @param print_result Whether to print the manipulation results. Default is FALSE. If set to TRUE, it will print the results of each manipulation step.
#' @param feature The features to be manipulated. If is_matrix is set to TRUE, feature will be ignored.
#' @param is_matrix Whether the data is in matrix format. Default is FALSE. If set to TRUE, the function will treat data as a matrix and extract the column names as features.
#'
#' @return filtered features
#' @export
#' @author Dongqiang Zeng
#' @examples

feature_manipulation<-function(data, feature, is_matrix = FALSE, print_result = FALSE){

  if(is_matrix){
    data<-as.data.frame(t(data))
    feature<-colnames(data)
  }

  data<-as.data.frame(data)
  #remove NA variables
  if(print_result)  print(paste(">>> Is NA exist: ",  sum(is.na(data[,feature]))))

  if (sum(is.na(data))>0) {
    nn<-as.data.frame(t(data[,feature]))
    delete_vars<- rownames(nn)[!complete.cases(nn)]
    feature<-feature[!feature%in%delete_vars]
  }
  #####################################

  #remove non-numeric variables
  if(print_result){
    print(paste0(">>>> Is nonnumeric variables exist ? >>>>"))
    print(summary(sapply(data[,feature], mode)!="numeric"))
  }
  fea_class<-as.vector(sapply(data[,feature], mode)=="numeric")
  feature<-feature[fea_class]

  #remove infinite variables

  if(print_result){
    print(paste0(">>>> Is -Inf variables exist ? >>>>"))
    print(summary(lapply(data[,feature],function(x) min(x))==-Inf))
  }

  fea_class<-as.vector(lapply(data[,feature],function(x) min(x))==-Inf)
  feature<-feature[!fea_class]

  if(print_result){
    print(paste0(">>>> Is Inf variables exist ? >>>>"))
    print(summary(lapply(data[,feature],function(x) max(x))== Inf))
  }
  fea_class<-as.vector(lapply(data[,feature],function(x) max(x))== Inf)
  feature<-feature[!fea_class]

  #remove variables with same number
  sd<-apply(data[,feature],2,function(x) sd(x)==0)

  if(print_result){
    print(paste0(">>> Variables with sd = 0 :  "))
    print(summary(sd))
  }

  feature<-feature[!sd]

  return(feature)
}
