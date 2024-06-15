



#' Merge expression sets by row names
#'
#' This function merges two or three expression sets by their row names, ensuring that 
#' no duplicated gene symbols remain after the merge. It is particularly useful in 
#' bioinformatics workflows where combining multiple datasets into a single expression 
#' set is required for downstream analysis.
#'
#' @param eset1 The first expression set to be merged.
#' @param eset2 The second expression set to be merged.
#' @param eset3 The optional third expression set to be merged; default is NULL.
#'
#' @return A merged expression set from the input expression sets, with duplicate gene symbols removed.
#' @export
#' 
#' @author Dongqiang Zeng
#' 
#' @examples
#' # Assuming `eset1`, `eset2`, and `eset3` are available in the environment
#' # All three expression sets must have the same row names (gene symbols)
#' merged_eset <- merge_eset(eset1, eset2, eset3)
merge_eset <- function(eset1, eset2, eset3 = NULL){

  eset <- merge(eset1, eset2, by = "row.names", all = FALSE)
  eset <-IOBR::remove_duplicate_genes(eset = eset, column_of_symbol = "Row.names")

  if(!is.null(eset3)){

    eset <- merge(eset, eset3, by = "row.names", all = FALSE)
    eset <-IOBR::remove_duplicate_genes(eset = eset, column_of_symbol = "Row.names")
  }

  return(eset)
}
