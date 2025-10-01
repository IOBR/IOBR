#' Merge Expression Sets by Row Names
#'
#' Merges two or three expression sets by row names (gene symbols), removing duplicates.
#'
#' @param eset1 First expression set.
#' @param eset2 Second expression set.
#' @param eset3 Optional third expression set. Default is NULL.
#'
#' @return Merged expression set with duplicates removed.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Assuming eset1, eset2, eset3 have same row names
#' merged_eset <- merge_eset(eset1, eset2, eset3)
merge_eset <- function(eset1, eset2, eset3 = NULL) {
  eset <- merge(eset1, eset2, by = "row.names", all = FALSE)
  eset <- IOBR::remove_duplicate_genes(eset = eset, column_of_symbol = "Row.names")

  if (!is.null(eset3)) {
    eset <- merge(eset, eset3, by = "row.names", all = FALSE)
    eset <- IOBR::remove_duplicate_genes(eset = eset, column_of_symbol = "Row.names")
  }

  return(eset)
}
