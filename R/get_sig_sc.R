#' get_sig_sc
#' @description Get gene signature from single-cell differential analysis
#' @param deg Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
#' @param cluster Name of the column in which the clusters are located
#' @param gene Name of the column in which the markers are located
#' @param n Number of selected top ranked markers
#' @param avg_log2FC Name of the column in which the average log2FC values are located
#'
#' @return A list containing top n gene markers of each cell types
#' @export
#'
#' @examples
#' data("deg", package = "IOBR")
#' get_sig_sc(deg, cluster = "cluster", gene = "gene", avg_log2FC = "avg_log2FC", n = 100)
get_sig_sc <- function(deg, cluster = "cluster", gene = "gene", avg_log2FC = "avg_log2FC", n = 100) {
  # cluster <- !!sym(cluster)
  # avg_log2FC <- !!sym(avg_log2FC)
  deg <- as.data.frame(deg)
  deg <- deg %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n, avg_log2FC)
  feas <- split(deg, deg[, cluster])
  feas <- lapply(feas, function(x) as.data.frame(x))
  feas <- lapply(feas, function(x) as.character(x[, gene]))
  feas <- lapply(feas, function(x) x[1:n])
  return(feas)
}
