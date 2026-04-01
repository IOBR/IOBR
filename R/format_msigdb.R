#' Format Input Signatures from MSigDB
#'
#' @description
#' Reads a GMT file from MSigDB and converts it into a named list of gene sets
#' suitable for IOBR functions.
#'
#' @param gmt Character string. Path to a GMT file.
#' @param ont Character string. Name of the signature/set column in the parsed
#'   GMT table. Default is `"term"`.
#' @param gene Character string. Name of the gene column in the parsed GMT
#'   table. Default is `"gene"`.
#'
#' @return Named list of character vectors, where each element contains the
#'   genes belonging to one signature.
#'
#' @export
#'
#' @examples
#' \donttest{
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
#' }
format_msigdb <- function(gmt, ont = "term", gene = "gene") {
  if (!file.exists(gmt)) {
    cli::cli_abort("File {.file {gmt}} does not exist")
  }

  sig <- clusterProfiler::read.gmt(gmt)
  sig <- as.data.frame(sig)

  colnames(sig)[colnames(sig) == ont] <- "ont"
  colnames(sig)[colnames(sig) == gene] <- "gene"

  sig_list <- split(sig[, c("ont", "gene")], sig$ont)
  sig_list <- lapply(sig_list, function(x) x$gene)

  sig_list
}
