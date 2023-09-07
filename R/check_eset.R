


#' check_eset - Check integrity and outliers of eset dataset
#'
#' @description This function is used to check the integrity and outliers of the given eset dataset, providing corresponding warnings and print results. It detects the presence of missing values, infinite values, and features with zero standard deviation.
#' @param eset The eset dataset to be checked.
#' @param print_result Whether to print the check results. Default is FALSE. If set to TRUE, it will print the results of each check.
#' @param estimate_sd Whether to estimate the standard deviation. Default is FALSE. If set to TRUE, it will estimate the standard deviation for each feature and provide warnings and print results accordingly.
#'
#' @return
#' @export
#'
#' @examples
check_eset<-function(eset, print_result = FALSE , estimate_sd = FALSE){



  if(print_result){
    message(paste("<<< Is NA exist ? >>> "))
    print(sum(is.na(eset)))

  }

  if(sum(is.na(eset))>0)  warning(paste0("There are some missing values in the matrix, which may affect the score calculation. You can set parameter 'adjust_eset' as TRUE to avoid these effects"))


  teset<-as.data.frame(t(eset))

  if(print_result){
    message(paste0("<<< Is -Inf feature exist ? >>>"))
    print(summary(lapply(teset,function(x) min(x))==-Inf))
  }


  if(print_result){
    message(paste0("<<< Is Inf feature exist ? >>>"))
    print(summary(lapply(teset,function(x) max(x))==Inf))
  }

  if(min(eset, na.rm = TRUE)==-Inf | max(eset)== Inf, na.rm = TRUE)  warning(paste0("There are some infinite values in the matrix, which may affect the score calculation. You can set parameter 'adjust_eset' as TRUE to avoid these effects"))


  if(estimate_sd){
    sd<-apply(eset,1,function(x) sd(x)==0)

    if(print_result){
      message(paste0("<<< Features have sd = 0 >>> "))
      print(summary(sd))
    }

    if(nlevels(as.factor(sd))>1) warning(paste0("Some vairables in the matrix have no variance between samples, which may affect the score calculation. You can set parameter 'adjust_eset' as TRUE to avoid these effects"))

  }
}
