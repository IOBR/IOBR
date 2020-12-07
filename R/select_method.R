


#' Choose which method to conduct following analysis
#'
#' @param data signature matrix
#' @param method choose method to analysis, contain three methodologies : PCA, ssGSEA, zScore
#'
#' @return filtered signature score matrix
#' @export
#'
#' @examples
select_method<-function(data, method = "ssGSEA" ){


  method<-tolower(method)

  if (method=="ssgsea") {
    data<-data[,-grep(colnames(data),pattern = "_PCA")]
    data<-data[,-grep(colnames(data),pattern = "_zscore")]
    colnames()<-gsub(colnames(data),pattern = "_ssGSEA",replacement = "")
  }
  if (method=="pca") {
    data<-data[,-grep(colnames(data),pattern = "_ssGSEA")]
    data<-data[,-grep(colnames(data),pattern = "_zscore")]
    colnames(data)<-gsub(colnames(data),pattern = "_PCA",replacement = "")
  }
  if (method=="zscore") {
    data<-data[,-grep(colnames(data),pattern = "_ssGSEA")]
    data<-data[,-grep(colnames(data),pattern = "_PCA")]
    colnames(data)<-gsub(colnames(data),pattern = "_zscore",replacement = "")
  }
  return(data)
}
#############################################
