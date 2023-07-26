




#' Scale matrix
#'
#' @param matrix matrix
#' @param log2matrix default is TRUE
#' @param manipulate
#'
#' @return
#' @export
#'
#' @examples
scale_matrix<-function(matrix, log2matrix = T, manipulate = FALSE){

  if(log2matrix){
    matrix<-IOBR::log2eset(matrix+1)
  }
  matrix<-as.data.frame(t(matrix))
  matrix<-scale(matrix,center = T,scale = T)
  matrix<-as.data.frame(t(matrix))

  if(manipulate){
    feas<-feature_manipulation(data = matrix, feature = rownames(matrix), is_matrix = T)
    matrix<- matrix[rownames(matrix)%in%feas, ]
  }
  return(matrix)
}
