#' Differential Expression Analysis
#'
#' @description
#' Performs differential expression analysis on gene expression data using either
#' DESeq2 or limma. Includes pre-processing steps like filtering low count data,
#' and calculates fold changes and adjusted p-values. Optionally generates
#' volcano plots and heatmaps.
#'
#' @param eset A matrix of gene expression data where rows represent genes and
#'   columns represent samples.
#' @param annotation Optional data frame for mapping gene IDs to gene names.
#'   Default is `NULL`.
#' @param pdata A data frame containing sample information and grouping labels.
#' @param group_id Character string specifying the column name in `pdata` containing
#'   grouping labels. Default is `"group"`.
#' @param pdata_id Character string specifying the column name in `pdata` for sample IDs.
#'   Default is `"ID"`.
#' @param array Logical indicating whether to perform quantile normalization.
#'   Default is `FALSE`.
#' @param method Character string specifying the method: `"DESeq2"` or `"limma"`.
#'   Default is `"DESeq2"`.
#' @param contrast Character vector of length 2 specifying contrast groups.
#'   Default is `c("High", "Low")`.
#' @param padj_cutoff Numeric cutoff for adjusted p-values. Default is `0.01`.
#' @param logfc_cutoff Numeric log2 fold change cutoff. Default is `0.5`.
#' @param volcano_plot Logical indicating whether to generate a volcano plot.
#'   Default is `FALSE`.
#' @param col_volcano Integer specifying color index for volcano plot. Default is `1`.
#' @param heatmap Logical indicating whether to generate a heatmap. Default is `TRUE`.
#' @param col_heatmap Integer specifying color index for heatmap. Default is `1`.
#' @param path Character string for output directory. Default is `NULL`.
#' @param parallel Logical indicating whether to run in parallel. Default is `FALSE`.
#' @param id_anno Character string specifying the identifier column in annotation.
#'   Default is `NULL`.
#'
#' @return Data frame containing differentially expressed genes with statistics
#'   including log2 fold changes and adjusted p-values.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' eset_stad <- load_data("eset_stad")
#' stad_group <- load_data("stad_group")
#' deg <- iobr_deg(
#'   eset = eset_stad, pdata = stad_group,
#'   group_id = "subtype", pdata_id = "ID", array = FALSE,
#'   method = "DESeq2", contrast = c("EBV", "GS"),
#'   path = file.path(tempdir(), "STAD")
#' )
#' head(deg)
#' }
iobr_deg <- function(eset,
                     annotation = NULL,
                     id_anno = NULL,
                     pdata,
                     group_id = "group",
                     pdata_id = "ID",
                     array = FALSE,
                     method = c("DESeq2", "limma"),
                     contrast = c("High", "Low"),
                     path = NULL,
                     padj_cutoff = 0.01,
                     logfc_cutoff = 0.5,
                     volcano_plot = FALSE,
                     col_volcano = 1,
                     heatmap = TRUE,
                     col_heatmap = 1,
                     parallel = FALSE) {
  method <- rlang::arg_match(method)

  # Input validation
  if (is.null(eset)) {
    cli::cli_abort("{.arg eset} cannot be NULL")
  }
  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort("{.arg eset} must be a matrix or data frame")
  }
  if (nrow(eset) == 0 || ncol(eset) == 0) {
    cli::cli_abort("{.arg eset} is empty")
  }
  if (is.null(pdata)) {
    cli::cli_abort("{.arg pdata} cannot be NULL")
  }
  pdata <- as.data.frame(pdata)
  if (!pdata_id %in% colnames(pdata)) {
    cli::cli_abort("Column {.val {pdata_id}} not found in pdata")
  }
  if (!group_id %in% colnames(pdata)) {
    cli::cli_abort("Column {.val {group_id}} not found in pdata")
  }
  if (length(contrast) != 2) {
    cli::cli_abort("{.arg contrast} must be a character vector of length 2")
  }

  # Setup output path
  abspath <- NULL
  if (!is.null(path)) {
    folder_info <- creat_folder(path)
    abspath <- folder_info$abspath
  }

  # Safely rename columns
  if (pdata_id != "ID") {
    colnames(pdata)[colnames(pdata) == pdata_id] <- "ID"
  }
  if (group_id != "deg_group") {
    colnames(pdata)[colnames(pdata) == group_id] <- "deg_group"
  }

  cli::cli_alert_info("Matching grouping information and expression matrix")

  if (group_id == "group3") {
    pdata <- pdata[pdata$deg_group != "Middle", ]
  }

  pdata <- pdata[!is.na(pdata$deg_group) & pdata$deg_group != "NA", ]
  pdata <- pdata[pdata$ID %in% colnames(eset), ]
  eset <- eset[, colnames(eset) %in% pdata$ID]
  pdata <- pdata[match(colnames(eset), pdata$ID), ]

  if (array) {
    rlang::check_installed("preprocessCore")
    eset <- preprocessCore::normalize.quantiles(as.matrix(eset), keep.names = TRUE)
  }

  if (method == "DESeq2") {
    DEG <- .run_deseq2(
      eset, pdata, contrast, annotation, id_anno,
      padj_cutoff, logfc_cutoff, parallel
    )
  } else {
    DEG <- .run_limma(eset, pdata, contrast, padj_cutoff, logfc_cutoff)
  }

  # Save results
  if (!is.null(abspath)) {
    save(DEG, file = paste0(abspath, "1-DEGs.RData"))
    csv_file <- paste0(abspath, "2-DEGs.csv")
    utils::write.csv(DEG, file = csv_file, row.names = FALSE)
    cli::cli_alert_success("DEG results written to: {.path {csv_file}}")
  }

  DEG
}

#' Run DESeq2 analysis
#' @keywords internal
#' @noRd
.run_deseq2 <- function(eset, pdata, contrast, annotation, id_anno,
                        padj_cutoff, logfc_cutoff, parallel) {
  rlang::check_installed("DESeq2")

  cli::cli_alert_info("Using DESeq2 for RNA-seq differential analysis")
  cli::cli_alert_warning("Ensure {.arg eset} is a count expression matrix")

  eset <- round(eset, 0)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = eset,
    colData = pdata,
    design = ~deg_group
  )

  # Filter low count genes
  dds <- dds[rowSums(DESeq2::counts(dds)) > ncol(eset) / 5, ]
  dds <- DESeq2::DESeq(dds, parallel = parallel)

  contrast <- c("deg_group", contrast)
  res_tidy <- DESeq2::results(dds, tidy = TRUE, contrast = contrast)
  res_tidy <- as.data.frame(res_tidy)

  cli::cli_alert_info("Differential analysis summary:")
  cli::cli_alert_info("  Adj.p < 0.001: {sum(res_tidy$padj < 0.001, na.rm = TRUE)}")
  cli::cli_alert_info("  Adj.p < 0.05: {sum(res_tidy$padj < 0.05, na.rm = TRUE)}")
  cli::cli_alert_info("  Adj.p < 0.1: {sum(res_tidy$padj < 0.1, na.rm = TRUE)}")
  cli::cli_alert_info("  Adj.p < 0.25: {sum(res_tidy$padj < 0.25, na.rm = TRUE)}")

  DEG <- res_tidy[order(res_tidy$padj, decreasing = FALSE), ]

  # Add annotation
  if (!is.null(annotation)) {
    colnames(annotation)[colnames(annotation) == id_anno] <- "id"
    DEG <- merge(DEG, annotation, by.x = "row", by.y = "id", all = FALSE)
  } else {
    cli::cli_alert_info("Using built-in anno_grch38 for annotation")
    anno_grch38 <- load_data("anno_grch38")
    DEG <- merge(DEG, anno_grch38, by.x = "row", by.y = "id", all = FALSE)
  }

  DEG <- .classify_degs(DEG, padj_cutoff, logfc_cutoff)

  # Get group IDs
  level1 <- as.character(contrast[2])
  level2 <- as.character(contrast[3])
  aa <- as.character(pdata[pdata$deg_group == level1, "ID"])
  aa <- aa[aa %in% colnames(eset)]
  bb <- as.character(pdata[pdata$deg_group == level2, "ID"])
  bb <- bb[bb %in% colnames(eset)]

  cli::cli_alert_info("Group 1 = {level1}: {length(aa)} samples")
  cli::cli_alert_info("Group 2 = {level2}: {length(bb)} samples")

  if (length(aa) <= 20 && length(bb) <= 20) {
    cli::cli_alert_info("Group 1 samples: {paste(aa, collapse = ', ')}")
    cli::cli_alert_info("Group 2 samples: {paste(bb, collapse = ', ')}")
  }

  # Calculate mean counts
  rlang::check_installed("preprocessCore")
  eset2 <- preprocessCore::normalize.quantiles(as.matrix(eset), keep.names = TRUE)
  eset2 <- tibble::rownames_to_column(as.data.frame(eset2), var = "ID")

  meancounts <- eset2 %>%
    dplyr::mutate(mean_group1 = rowSums(.[, aa, drop = FALSE]) / length(aa)) %>%
    dplyr::mutate(mean_group2 = rowSums(.[, bb, drop = FALSE]) / length(bb)) %>%
    dplyr::select("ID", "mean_group1", "mean_group2")

  DEG <- merge(DEG, meancounts, by.x = "row", by.y = "ID", all = FALSE)
  colnames(DEG)[colnames(DEG) == "mean_group1"] <- level1
  colnames(DEG)[colnames(DEG) == "mean_group2"] <- level2

  DEG <- tibble::as_tibble(DEG)
  DEG[order(DEG$padj, decreasing = FALSE), ]
}

#' Run limma analysis
#' @keywords internal
#' @noRd
.run_limma <- function(eset, pdata, contrast, padj_cutoff, logfc_cutoff) {
  rlang::check_installed("limma")

  cli::cli_alert_info("Using limma for array differential analysis")

  contrast_full <- c("deg_group", contrast)
  keep_groups <- contrast[1:2]
  pdata <- pdata[pdata$deg_group %in% keep_groups, , drop = FALSE]
  pdata <- pdata[!is.na(pdata$deg_group), , drop = FALSE]
  eset <- eset[, colnames(eset) %in% pdata$ID, drop = FALSE]
  pdata <- pdata[match(colnames(eset), pdata$ID), , drop = FALSE]

  if (!all(keep_groups %in% unique(as.character(pdata$deg_group)))) {
    available <- paste(unique(as.character(pdata$deg_group)), collapse = ", ")
    cli::cli_abort("Selected groups not found. Available: {available}")
  }

  pdata$deg_group <- ifelse(pdata$deg_group == contrast[1], "group1", "group2")

  cli::cli_alert_info("Group 1 = {contrast[1]}")
  cli::cli_alert_info("Group 2 = {contrast[2]}")

  design <- stats::model.matrix(~ 0 + factor(pdata$deg_group))
  colnames(design) <- levels(factor(pdata$deg_group))
  rownames(design) <- colnames(eset)

  contrast.matrix <- limma::makeContrasts(group1 - group2, levels = design)
  fit <- limma::lmFit(eset, design)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  DEG <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

  DEG <- tibble::rownames_to_column(DEG, var = "symbol")
  colnames(DEG)[colnames(DEG) == "logFC"] <- "log2FoldChange"
  colnames(DEG)[colnames(DEG) == "adj.P.Val"] <- "padj"
  colnames(DEG)[colnames(DEG) == "P.Value"] <- "pvalue"

  DEG <- .classify_degs(DEG, padj_cutoff, logfc_cutoff)

  # Get mean counts
  aa <- as.character(pdata[pdata$deg_group == "group1", "ID"])
  aa <- aa[aa %in% colnames(eset)]
  bb <- as.character(pdata[pdata$deg_group == "group2", "ID"])
  bb <- bb[bb %in% colnames(eset)]

  eset2 <- tibble::rownames_to_column(as.data.frame(eset), var = "ID")
  meancounts <- eset2 %>%
    dplyr::mutate(mean_group1 = rowSums(.[, aa, drop = FALSE]) / length(aa)) %>%
    dplyr::mutate(mean_group2 = rowSums(.[, bb, drop = FALSE]) / length(bb)) %>%
    dplyr::select("ID", "mean_group1", "mean_group2")

  DEG <- merge(DEG, meancounts, by.x = "symbol", by.y = "ID", all = FALSE)
  colnames(DEG)[colnames(DEG) == "mean_group1"] <- contrast[1]
  colnames(DEG)[colnames(DEG) == "mean_group2"] <- contrast[2]

  tibble::as_tibble(DEG[order(DEG$padj, decreasing = FALSE), ])
}

#' Classify DEGs
#' @keywords internal
#' @noRd
.classify_degs <- function(DEG, padj_cutoff, logfc_cutoff) {
  DEG$sigORnot <- ifelse(
    DEG$log2FoldChange > logfc_cutoff & DEG$padj < padj_cutoff, "Up_regulated",
    ifelse(DEG$log2FoldChange < -logfc_cutoff & DEG$padj < padj_cutoff, "Down_regulated", "NOT")
  )
  DEG$label <- ifelse(
    abs(DEG$log2FoldChange) > logfc_cutoff & DEG$padj < padj_cutoff, "Both",
    ifelse(DEG$padj < padj_cutoff, "Significant",
      ifelse(abs(DEG$log2FoldChange) >= logfc_cutoff, paste0("log2FC >= ", logfc_cutoff), "NOT")
    )
  )
  DEG
}
