#' Load IOBR Datasets
#'
#' @param name A dataset name.
#'
#' @returns Dataset, typically `list` or `data.frame`.
#' @export
#'
#' @examples
#' load_data("deg")
load_data <- function(name) {
  stopifnot(length(name) == 1)
  sysdata_names <- c(
    "cancer_type_genes", "cellmarkers", "common_genes", "go_bp",
    "go_cc", "go_mf", "hallmark", "immuneCuratedData", "imvigor210_eset",
    "ips_gene_set", "kegg", "length_ensembl", "mcp_genes", "mcp_probesets",
    "melanoma_data", "mRNA_cell_default", "msig_immune", "msig_sc",
    "mus_human_gene_symbol", "onco_sig", "palette1", "palette2",
    "palette3", "palette4", "panel_for_gene", "panel_for_signature",
    "patterns_to_na", "pdata_acrg", "pdata_GSE63557", "pdata_sig_tme",
    "pdata_sig_tme_binary", "pdata_tme_binary", "PurityDataAffy",
    "quantiseq_data", "reactome", "SI_geneset", "sig_excel", "signature_collection_citation",
    "signature_metabolism", "signature_sc", "signature_tme", "signature_tumor",
    "tcga_stad_var", "xCell.data"
  )

  data_names <- utils::data(package = "IOBR")$results[, "Item"]

  if (name %in% sysdata_names) {
    eval(parse(text = name))
  } else if (name %in% data_names) {
    data(list = name, package = "IOBR", envir = environment())
    get(name)
  } else {
    print("Available datasets:")
    print(c(sysdata_names, data_names))
    stop("Input dataset name not found, availables see above")
  }
}
