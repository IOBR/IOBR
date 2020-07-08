


#' Choose which method to conduct analysis
#'
#' @param score TME signature score matrix
#' @param choosemethod choosemethod to analysis, contain three methodologies : PCA, ssGSEA, zScore
#'
#' @return filterd TME signature score matrix
#' @export
#'
#' @examples
select_method<-function(score = score, choosemethod = "ssGSEA" ){

  if (choosemethod=="ssGSEA") {
    score<-score[,-grep(colnames(score),pattern = "_PCA")];
    score<-score[,-grep(colnames(score),pattern = "_zscore")]
    colnames(score)<-gsub(colnames(score),pattern = "_ssGSEA",replacement = "")
  }
  if (choosemethod=="PCA") {
    score<-score[,-grep(colnames(score),pattern = "_ssGSEA")];
    score<-score[,-grep(colnames(score),pattern = "_zscore")]
    colnames(score)<-gsub(colnames(score),pattern = "_PCA",replacement = "")
  }
  if (choosemethod=="zScore") {
    score<-score[,-grep(colnames(score),pattern = "_ssGSEA")];
    score<-score[,-grep(colnames(score),pattern = "_PCA")]
    colnames(score)<-gsub(colnames(score),pattern = "_zscore",replacement = "")
  }
  return(score)
}
#############################################
