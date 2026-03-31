#' Extract Data Frame from Seurat Object
#'
#' @description
#' Extracts and combines a data frame with cells as rows and features as columns
#' from Seurat assay data. Supports multiple assays and optional metadata
#' integration.
#'
#' @param sce Seurat object.
#' @param vars Character vector of feature names to extract. If `NULL`, all
#'   features are extracted.
#' @param assay Character vector specifying assay(s) to pull data from.
#' @param slot Character string specifying the assay data slot. Default is
#'   `"scale.data"`.
#' @param combine_meta_data Logical indicating whether to combine metadata
#'   with the extracted data frame. Default is `TRUE`.
#'
#' @return Data frame with cells as rows and features as columns.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \dontrun{
#' pbmc_small <- load_data("pbmc_small")
#' vars <- c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3")
#' eset <- extract_sc_data(sce = pbmc_small, vars = vars, assay = "RNA")
#' }
extract_sc_data <- function(sce, vars = NULL, assay, slot = "scale.data", combine_meta_data = TRUE) {
  rlang::check_installed(c("Seurat", "SeuratObject"))

  if (!inherits(sce, "Seurat")) {
    cli::cli_abort("{.arg sce} must be a Seurat object")
  }

  exist <- Seurat::Assays(sce)
  cli::cli_alert_info("Available assays: {.val {exist}}")

  assay <- assay[assay %in% exist]
  if (length(assay) == 0) {
    cli::cli_abort("No valid assay found in object")
  }

  eset_cbind <- data.frame(
    ID = rownames(sce@meta.data),
    index = seq_len(nrow(sce@meta.data))
  )

  for (i in seq_along(assay)) {
    method <- assay[i]
    Seurat::DefaultAssay(sce) <- method

    eset <- SeuratObject::GetAssayData(sce, assay = method, slot = slot)

    if (!is.null(vars)) {
      feas <- rownames(eset)[rownames(eset) %in% unique(vars)]
      if (length(feas) == 0) {
        cli::cli_abort("Required variables not found in expression matrix")
      }
      eset <- eset[feas, , drop = FALSE]
    }

    eset <- as.data.frame(t(as.matrix(eset)))

    if (!is.null(vars) && length(vars) == 1) {
      eset <- as.data.frame(eset)
      rownames(eset) <- vars
      eset <- t(eset)
      eset <- as.data.frame(eset)
    }

    eset <- tibble::rownames_to_column(eset, var = "ID")

    if (length(assay) > 1) {
      colnames(eset)[2:ncol(eset)] <- paste0(
        colnames(eset)[2:ncol(eset)], "_", method
      )
    }

    eset <- as.data.frame(eset)

    if (length(assay) == 1) {
      eset_cbind <- eset
    } else {
      eset_cbind <- dplyr::inner_join(eset_cbind, eset, by = "ID")
    }
  }

  if (combine_meta_data) {
    cli::cli_alert_info("Merging metadata...")
    meta.data <- tibble::rownames_to_column(sce@meta.data, var = "ID")
    eset_cbind <- dplyr::inner_join(meta.data, eset_cbind, by = "ID")
  }

  eset_cbind
}
