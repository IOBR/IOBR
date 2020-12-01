


#' Format input signatures from MsgiDB
#'
#' @param sig `signature` data frame downloaded from `https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H`
#' @param ont column name of signatures
#' @param gene column name of gene
#'
#' @return signature list
#' @export
#'
#' @examples
format_msigdb<-function(gmt, ont = "ont", gene = "gene"){

  sig <- clusterProfiler::read.gmt(paste0(gmt))
  ######################################
  sig<-as.data.frame(sig)
  colnames(sig)[which(colnames(sig)==ont)]<-"ont"
  colnames(sig)[which(colnames(sig)==gene)]<-"gene"
  sig_list <- select(sig, ont, gene) %>%
    as.data.frame %>%
    split(., .$ont) %>%
    lapply(., function(x)(x$gene))
  return(sig_list)
}
