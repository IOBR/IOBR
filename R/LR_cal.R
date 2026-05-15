#' Calculate Ligand-Receptor Interaction Scores
#'
#' @description
#' Quantifies ligand-receptor interactions in the tumor microenvironment from
#' bulk gene expression data using the easier package. This function processes
#' raw counts or TPM data and computes interaction scores for each sample.
#'
#' @param eset Gene expression matrix with genes as rows and samples as columns.
#' @param data_type Type of input data. Options are `"count"` or `"tpm"`.
#'   If `"count"`, data will be converted to TPM before analysis.
#' @param id_type Type of gene identifier. Default is `"ensembl"`.
#' @param cancer_type Character string specifying the cancer type for easier.
#'   Default is `"pancan"` for pan-cancer analysis.
#'
#' @return Data frame containing ligand-receptor interaction scores with
#'   sample IDs as row names.
#'
#' @references
#' Lapuente-Santana, van Genderen, M., Hilbers, P., Finotello, F., & Eduati, F. (2021).
#' Interpretable systems biomarkers predict response to immune-checkpoint inhibitors.
#' Patterns (New York, N.Y.), 2(8), 100293. https://doi.org/10.1016/j.patter.2021.100293
#'
#' @export
#'
#' @examples
#' # LR_cal requires HGNC gene symbols as rownames
#' # Create a simple example with gene symbols
#' example_genes <- c(
#'   "TGFB1", "EGFR", "VEGFA", "PDGFB", "FGF2", "CXCL12",
#'   "CXCR4", "IL6", "IL6R", "TNF", "TNFRSF1A", "IFNG"
#' )
#' sim_eset <- as.data.frame(matrix(
#'   rnorm(length(example_genes) * 10, mean = 5, sd = 2),
#'   nrow = length(example_genes), ncol = 10
#' ))
#' rownames(sim_eset) <- example_genes
#' colnames(sim_eset) <- paste0("Sample", 1:10)
#' \dontrun{
#' if (requireNamespace("easier", quietly = TRUE)) {
#'   tryCatch({
#'     lr <- LR_cal(eset = sim_eset, data_type = "tpm")
#'     head(lr)
#'   }, error = function(e) {
#'     message("Example skipped: could not download ExperimentHub data")
#'   })
#' }
#' }
LR_cal <- function(eset, data_type = c("count", "tpm"), id_type = "ensembl", cancer_type = "pancan") {
  rlang::check_installed("easier")
  data_type <- rlang::arg_match(data_type)

  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort("{.arg eset} must be a matrix or data frame")
  }
  if (nrow(eset) == 0 || ncol(eset) == 0) {
    cli::cli_abort("{.arg eset} must have at least one row and one column")
  }

  if (data_type == "count") {
    cli::cli_alert_info("Converting count data to TPM...")
    eset <- count2tpm(countMat = eset, idType = id_type, source = "local")
  }

  eset <- as.matrix(eset)
  feas <- feature_manipulation(data = eset, feature = rownames(eset), is_matrix = TRUE)
  eset <- eset[rownames(eset) %in% feas, , drop = FALSE]

  if (nrow(eset) == 0) {
    cli::cli_abort("No valid features remaining after filtering")
  }

  res <- easier::compute_LR_pairs(RNA_tpm = eset, cancer_type = cancer_type, verbose = TRUE)

  tibble::rownames_to_column(res, var = "ID")
}
