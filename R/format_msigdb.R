#' Format input signatures from MSigDB
#'
#' Read a GMT file from MSigDB and convert it into a named list of gene sets.
#'
#' @param gmt Character string. Path to a GMT file.
#' @param ont Character string. Name of the signature/set column in the parsed GMT table.
#'   Default is \code{"term"}.
#' @param gene Character string. Name of the gene column in the parsed GMT table.
#'   Default is \code{"gene"}.
#'
#' @return A named list of character vectors, where each element contains the genes
#'   belonging to one signature.
#' @export
#'
#' @examples
#' # Create a minimal GMT file for the example
#' tf <- tempfile(fileext = ".gmt")
#' writeLines(
#'   c(
#'     "HALLMARK_TNFA_SIGNALING_VIA_NFKB\tNA\tTNF\tNFKB1\tNFKB2",
#'     "HALLMARK_P53_PATHWAY\tNA\tTP53\tMDM2\tCDKN1A"
#'   ),
#'   con = tf
#' )
#'
#' sig_list <- format_msigdb(tf, ont = "term", gene = "gene")
#' names(sig_list)
#' sig_list[[1]]
format_msigdb <- function(gmt, ont = "term", gene = "gene") {
  sig <- clusterProfiler::read.gmt(gmt)
  sig <- as.data.frame(sig)

  colnames(sig)[colnames(sig) == ont] <- "ont"
  colnames(sig)[colnames(sig) == gene] <- "gene"

  sig_list <- sig[, c("ont", "gene")] %>%
    split(.$ont) %>%
    lapply(function(x) x$gene)

  return(sig_list)
}
