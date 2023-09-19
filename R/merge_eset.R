



#' Merge expression sets by row names
#'
#' @param eset1 eset1
#' @param eset2 eset2
#' @param eset3 eset3, default is NULL
#'
#' @return
#' @export
#'
#' @examples
#' @author Dongqiang Zeng
merge_eset <- function(eset1, eset2, eset3 = NULL){
  
  eset <- merge(eset1, eset2, by = "row.names", all = FALSE)
  eset <-IOBR::remove_duplicate_genes(eset = eset, column_of_symbol = "Row.names")
  
  if(!is.null(eset3)){
    
    eset <- merge(eset, eset3, by = "row.names", all = FALSE)
    eset <-IOBR::remove_duplicate_genes(eset = eset, column_of_symbol = "Row.names")
  }
  
}