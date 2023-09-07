




#' scale_matrix - Scale and manipulate a matrix
#'
#' @description This function scales a given matrix and optionally performs additional manipulation. It can apply a logarithmic transformation (base 2) to the matrix, center and scale the values, and optionally manipulate the features based on a specified criteria.
#' @param matrix The matrix to be scaled and manipulated.
#' @param log2matrix Whether to apply a logarithmic transformation to the matrix. Default is TRUE. If set to TRUE, the matrix will be transformed using the log2eset function from the IOBR package.
#' @param manipulate Whether to perform additional manipulation on the matrix. Default is FALSE. If set to TRUE, the function will use the feature_manipulation function to manipulate the matrix based on specified criteria.
#'
#' @return
#' @export
#'
#' @author Dongqiang Zeng
#' @examples
scale_matrix<-function(matrix, log2matrix = TRUE, manipulate = TRUE){

  if(log2matrix){
    matrix<- IOBR::log2eset(matrix+1)
  }
  matrix<-as.data.frame(t(matrix))
  matrix<-scale(matrix,center = T,scale = T)
  matrix<-as.data.frame(t(matrix))

  if(manipulate){
    feas<-feature_manipulation(data = matrix, feature = rownames(matrix), is_matrix = TRUE, print_result = TRUE)
    matrix<- matrix[rownames(matrix)%in%feas, ]
  }
  return(matrix)
}
