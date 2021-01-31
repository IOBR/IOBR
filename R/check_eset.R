


#' Check the value of gene expression set
#'
#' @param eset gene expression set
#' @param print_result logical, default is FALSE
#'
#' @return
#' @export
#'
#' @examples
check_eset<-function(eset, print_result = F){
  
  
  
  if(print_result){
    message(paste("<<< Is NA exist ? >>> "))
    print(sum(is.na(eset)))
    
  }  
  
  if(sum(is.na(eset))>0)  warning(paste0("There are some missing values in the gene expression matrix, which may affect the score calculation. You can set parameter 'fea_engine' as TRUE to avoid these effects"))
  
  
  teset<-as.data.frame(t(eset))
  
  if(print_result){
    message(paste0("<<< Is -Inf feature exist ? >>>"))
    print(summary(lapply(teset,function(x) min(x))==-Inf))
  }
  
  
  if(print_result){
    message(paste0("<<< Is Inf feature exist ? >>>"))
    print(summary(lapply(teset,function(x) max(x))==Inf))
  }
  
  if(min(eset)==-Inf | max(eset)== Inf)  warning(paste0("There are some infinite values in the gene expression matrix, which may affect the score calculation. You can set parameter 'fea_engine' as TRUE to avoid these effects"))
  
  
  sd<-apply(eset,1,function(x) sd(x)==0)
  
  if(print_result){
    message(paste0("<<< Features have sd = 0 >>> "))
    print(summary(sd))
  }
  
  if(nlevels(as.factor(sd))>1) warning(paste0("Some vairables in the gene expression matrix have no variance between samples, which may affect the score calculation. You can set parameter 'fea_engine' as TRUE to avoid these effects"))
  
  
}
