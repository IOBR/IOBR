#' Extract Top Marker Genes from Single-Cell Differential Results
#'
#' @description
#' Selects the top N marker genes per cluster from a ranked differential
#' expression result table.
#'
#' @param deg Data frame or matrix. Ranked marker statistics.
#' @param cluster Character. Column name containing cluster identifiers.
#'   Default is `"cluster"`.
#' @param gene Character. Column name containing gene identifiers.
#'   Default is `"gene"`.
#' @param n Integer. Number of top markers per cluster. Default is 100.
#' @param avg_log2FC Character. Column name for average log2 fold change.
#'   Default is `"avg_log2FC"`.
#'
#' @return List of character vectors; each element contains the top N genes
#'   for a cluster.
#'
#' @export
#'
#' @examples
#' deg <- load_data("deg")
#' get_sig_sc(deg, cluster = "cluster", gene = "gene", avg_log2FC = "avg_log2FC", n = 100)
get_sig_sc <- function(deg, cluster = "cluster", gene = "gene", avg_log2FC = "avg_log2FC", n = 100) {
  if (!is.data.frame(deg) && !is.matrix(deg)) {
    cli::cli_abort("{.arg deg} must be a data frame or matrix")
  }
  if (!cluster %in% colnames(deg)) {
    cli::cli_abort("Column {.val {cluster}} not found in {.arg deg}")
  }
  if (!gene %in% colnames(deg)) {
    cli::cli_abort("Column {.val {gene}} not found in {.arg deg}")
  }
  if (!avg_log2FC %in% colnames(deg)) {
    cli::cli_abort("Column {.val {avg_log2FC}} not found in {.arg deg}")
  }

  deg <- as.data.frame(deg)

  deg_top <- deg %>%
    dplyr::group_by(.data[[cluster]]) %>%
    dplyr::top_n(n, .data[[avg_log2FC]])

  feas <- split(deg_top, deg_top[[cluster]])
  feas <- lapply(feas, function(x) {
    genes <- as.character(x[[gene]])
    head(genes, n)
  })

  feas
}
