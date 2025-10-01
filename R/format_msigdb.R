#' Format input signatures from MsgiDB
#'
#' @param gmt `signature` data frame downloaded from `https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H`
#' @param ont column name of signatures
#' @param gene column name of gene
#'
#' @return signature list
#' @export
#'
#' @examples
#' # download gmt through: https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt
#' format_msigdb("E:/0-0-0-My Pipeline/00-Signature-collection/h.all.v2023.2.Hs.symbols.gmt", ont = "term")
#'
format_msigdb <- function(gmt, ont = "term", gene = "gene") {
  sig <- clusterProfiler::read.gmt(paste0(gmt))
  ######################################
  sig <- as.data.frame(sig)
  colnames(sig)[which(colnames(sig) == ont)] <- "ont"
  colnames(sig)[which(colnames(sig) == gene)] <- "gene"
  sig_list <- select(sig, ont, gene) %>%
    as.data.frame() %>%
    split(., .$ont) %>%
    lapply(., function(x) (x$gene))
  return(sig_list)
}
