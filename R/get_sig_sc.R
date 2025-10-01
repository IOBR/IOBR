#' Extract Top Marker Genes from Single-Cell Differential Results
#'
#' Selects the top N marker genes per cluster from a ranked differential expression result table.
#'
#' @param deg Data frame or matrix. Ranked marker statistics (e.g., p-value, log2FC, etc.).
#' @param cluster Character. Column name containing cluster identifiers. Default "cluster".
#' @param gene Character. Column name containing gene identifiers. Default "gene".
#' @param n Integer. Number of top markers per cluster. Default 100.
#' @param avg_log2FC Character. Column name for average log2 fold change. Default "avg_log2FC".
#'
#' @return List of character vectors; each element contains the top N genes for a cluster.
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
