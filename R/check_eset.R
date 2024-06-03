


#' check_eset - Check integrity and outliers of expression set
#'
#' @description This function is used to check the integrity and outliers of the given expression set, providing corresponding warnings and print results. It detects the presence of missing values, infinite values, and features with zero standard deviation.
#' @param eset The matrix to be checked.
#' @param print_result Whether to print the check results. Default is FALSE. If set to TRUE, it will print the results of each check.
#' @param estimate_sd Whether to estimate the standard deviation. Default is FALSE. If set to TRUE, it will estimate the standard deviation for each feature and provide warnings and print results accordingly.
#'
#' @return NULL. The function is called for its side effects of checking and optionally printing warnings about the integrity of the expression set.
#' @export
#'
#' @author Dongqiang Zeng
#' @examples
#' # Loading TCGA-STAD expresion data(raw count matrix)
#' data("eset_stad", package = "IOBR")
#' # transform count data to tpm
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' check_eset(eset)
check_eset<-function(eset, print_result = FALSE , estimate_sd = FALSE){



  if(print_result){
    message(paste("<<< Do NA values exist ? >>> "))
    print(sum(is.na(eset)))

  }

  if(sum(is.na(eset))>0)  warning(paste0("There are some missing values in the matrix, which may affect the score calculation. You can set parameter 'adjust_eset' as TRUE to avoid these effects"))


  teset<-as.data.frame(t(eset))

  if(print_result){
    message(paste0("<<< Do -Inf features exist ? >>>"))
    print(summary(lapply(teset,function(x) min(x))==-Inf))
  }


  if(print_result){
    message(paste0("<<< Do Inf features exist ? >>>"))
    print(summary(lapply(teset,function(x) max(x))==Inf))
  }

  if(min(eset, na.rm = TRUE)== -Inf | max(eset, na.rm = TRUE)== Inf)  warning(paste0("There are infinite values in the matrix, which may affect the score calculation. You can set the 'adjust_eset' parameter as TRUE to avoid these effects."))


  if(estimate_sd){
    sd<-apply(eset,1,function(x) sd(x)==0)

    if(print_result){
      message(paste0("<<< Features have sd = 0 >>> "))
      print(summary(sd))
    }

    if(nlevels(as.factor(sd))>1) warning(paste0("Some variables in the matrix have no variance between samples, which may affect the score calculation. You can set the 'adjust_eset' parameter as TRUE to avoid these effects."))

  }
}
