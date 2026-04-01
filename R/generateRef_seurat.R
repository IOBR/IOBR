#' Generate Reference Matrix from Seurat Object
#'
#' @description
#' Generates reference gene expression data from a Seurat object by
#' identifying marker genes for each cell type and aggregating expression
#' data.
#'
#' @param sce Seurat object containing single-cell RNA-seq data.
#' @param celltype Character. Cell type column name in metadata.
#'   Default is `NULL` (uses default identity).
#' @param proportion Numeric. Proportion of cells to randomly select for
#'   analysis. Default is `NULL` (use all cells).
#' @param assay_deg Character. Assay for finding markers. Default is `"RNA"`.
#' @param slot_deg Character. Slot for finding markers. Default is `"data"`.
#' @param adjust_assay Logical. Whether to adjust assay for SCT. Default is
#'   `FALSE`.
#' @param assay_out Character. Assay for output. Default is `"RNA"`.
#' @param slot_out Character. Slot for output. Default is `"data"`.
#' @param verbose Logical. Print verbose messages. Default is `FALSE`.
#' @param only.pos Logical. Return only positive markers. Default is `TRUE`.
#' @param n_ref_genes Integer. Number of reference genes per cell type.
#'   Default is 50.
#' @param logfc.threshold Numeric. Log fold change threshold. Default is 0.15.
#' @param test.use Character. Statistical test for marker identification.
#'   Default is `"wilcox"`.
#'
#' @return Matrix containing aggregated expression data for reference genes.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("SeuratObject", quietly = TRUE)) {
#'   pbmc <- SeuratObject::pbmc_small
#'   sm <- generateRef_seurat(sce = pbmc, celltype = "groups", slot_out = "data")
#' }
#' }
generateRef_seurat <- function(sce, celltype = NULL, proportion = NULL,
                               assay_deg = "RNA", slot_deg = "data",
                               adjust_assay = FALSE, assay_out = "RNA",
                               slot_out = "data", verbose = FALSE,
                               only.pos = TRUE, n_ref_genes = 50,
                               logfc.threshold = 0.15, test.use = "wilcox") {
  rlang::check_installed("Seurat")

  if (!inherits(sce, "Seurat")) {
    cli::cli_abort("{.arg sce} must be a Seurat object")
  }

  cli::cli_alert_info("Assay used to find markers: {assay_deg}")

  if (!is.null(assay_deg)) {
    Seurat::DefaultAssay(sce) <- assay_deg
  }

  if (!is.null(celltype)) {
    cli::cli_alert_info("Idents of Seurat object is: {celltype}")
    Seurat::Idents(sce) <- celltype
    print(table(as.character(Seurat::Idents(sce))))
  } else {
    cli::cli_alert_info("Idents of Seurat object is: Default...")
    print(table(as.character(Seurat::Idents(sce))))
  }

  if (tolower(assay_deg) == "sct" && adjust_assay) {
    sce <- Seurat::PrepSCTFindMarkers(sce)
  }

  if (!is.null(proportion)) {
    input <- random_strata_cells(
      input = sce@meta.data,
      group = celltype,
      propotion = proportion
    )
    cell_input <- rownames(input)

    cli::cli_alert_info(
      "Subsetting cells in each cluster randomly: {proportion * 100}%..."
    )
    cli::cli_alert_info("Final training data:")
    print(table(as.factor(input[, celltype])))

    sce <- subset(sce, cells = as.character(cell_input))
  }

  cli::cli_alert_info("Find markers of each celltype...")

  sce.markers <- Seurat::FindAllMarkers(
    object = sce,
    slot = slot_deg,
    assay = assay_deg,
    features = NULL,
    only.pos = only.pos,
    min.pct = 0.25,
    thresh.use = 0.25,
    verbose = verbose,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    fc.name = "avg_log2FC"
  )

  print(utils::head(sce.markers))

  refgene <- sce.markers |>
    dplyr::group_by(.data$cluster) |>
    dplyr::top_n(n_ref_genes, .data$avg_log2FC)

  print(refgene)

  cli::cli_alert_info("Aggregating scRNAseq data...")

  if (is.null(slot_out)) slot_out <- "counts"

  cli::cli_alert_info(
    "orig.ident was set as group. User can define through parameter celltype..."
  )

  bulk <- Seurat::AggregateExpression(
    object = sce,
    assays = assay_out,
    slot = slot_out,
    group.by = celltype
  )
  bulk <- bulk[[1]]

  bulk <- bulk[rownames(bulk) %in% refgene$gene, , drop = FALSE]

  bulk
}
