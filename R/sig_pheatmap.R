#' Generate Heatmap for Signature Data
#'
#' @description
#' Creates a heatmap from signature data with grouping variables, offering
#' flexible options for colors, clustering, and output formats using
#' ComplexHeatmap.
#'
#' @param input Data frame with variables in columns.
#' @param feas Character vector. Feature names (columns) to include in heatmap.
#' @param group Character string. Column name for primary grouping variable.
#' @param group2 Character string or `NULL`. Optional secondary grouping variable.
#' @param group3 Character string or `NULL`. Optional tertiary grouping variable.
#' @param ID Character string. Column name for sample identifiers. Default is `"ID"`.
#' @param path Character string or `NULL`. Directory to save output files.
#'   Default creates `"Marker-heatmap-average"`.
#' @param cols1 Character vector or `"random"` or `"normal"`. Colors for primary group.
#'   Default is `"random"`.
#' @param seed Integer. Random seed for color generation. Default is `54321`.
#' @param show_col Logical indicating whether to display colors. Default is `FALSE`.
#' @param cluster_cols Logical indicating whether to cluster columns. Default is `TRUE`.
#' @param palette_for_heatmape Integer. Palette number for heatmap. Default is `6`.
#' @param scale.matrix Logical indicating whether to scale the matrix. Default is `TRUE`.
#' @param cellwidth Numeric. Width of each cell in points. Default is `1`.
#' @param cellheight Numeric. Height of each cell in points. Default is `9`.
#' @param fig.type Character string. File format for saving. Default is `"pdf"`.
#' @param width Numeric. Width of saved figure in inches. Default is `6`.
#' @param height Numeric or `NULL`. Height of saved figure in inches. Calculated if `NULL`.
#' @param file_name_prefix Character or numeric. Prefix for saved file name. Default is `1`.
#' @param cols2 Character vector or `"random"` or `"normal"`. Colors for secondary group.
#'   Default is `"random"`.
#' @param cols3 Character vector or `"random"` or `"normal"`. Colors for tertiary group.
#'   Default is `"random"`.
#' @param palette1 Integer. Palette for primary group. Default is `1`.
#' @param palette2 Integer. Palette for secondary group. Default is `2`.
#' @param palette3 Integer. Palette for tertiary group. Default is `3`.
#' @param show_colnames Logical indicating whether to show column names. Default is `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{p_anno}{Annotation data frame}
#'   \item{p_cols}{List of cluster colors}
#'   \item{plot}{ComplexHeatmap object}
#'   \item{eset}{Transformed expression matrix}
#' }
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' tcga_stad_sig <- load_data("tcga_stad_sig")
#' tcga_stad_pdata <- load_data("tcga_stad_pdata")
#' input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
#' feas <- grep("MCPcounter", colnames(input), value = TRUE)
#' sig_pheatmap(input = input, feas = feas, group = "subtype", scale.matrix = TRUE)
#' }
sig_pheatmap <- function(input,
                         feas,
                         group,
                         group2 = NULL,
                         group3 = NULL,
                         ID = "ID",
                         path = NULL,
                         cols1 = "random",
                         cols2 = "random",
                         cols3 = "random",
                         seed = 54321,
                         show_col = FALSE,
                         palette1 = 1,
                         palette2 = 2,
                         palette3 = 3,
                         cluster_cols = TRUE,
                         palette_for_heatmape = 6,
                         scale.matrix = TRUE,
                         cellwidth = 1,
                         cellheight = 9,
                         show_colnames = FALSE,
                         fig.type = "pdf",
                         width = 6,
                         height = NULL,
                         file_name_prefix = 1) {
  rlang::check_installed("ComplexHeatmap")
  rlang::check_installed("grid")

  # Create output directory
  file_store <- path %||% "Marker-heatmap-average"
  path_obj <- creat_folder(file_store)

  # Prepare input data
  input <- as.data.frame(input)
  feas <- feas[feas %in% colnames(input)]

  if (length(feas) == 0) {
    cli::cli_abort("No valid features found in input")
  }

  if (!ID %in% colnames(input)) {
    cli::cli_abort("ID column {.val {ID}} not found in input")
  }

  if (!group %in% colnames(input)) {
    cli::cli_abort("Group column {.val {group}} not found in input")
  }

  # Create expression matrix
  eset <- input[, c(ID, feas), drop = FALSE]
  rownames(eset) <- NULL
  eset <- tibble::column_to_rownames(eset, var = ID)

  if (scale.matrix) {
    eset <- scale(eset, scale = TRUE, center = TRUE)
  }
  eset <- t(eset)

  # Create annotation data frame
  anno_cols <- c(ID, group)
  if (!is.null(group2)) anno_cols <- c(anno_cols, group2)
  if (!is.null(group3)) anno_cols <- c(anno_cols, group3)

  anno <- input[, anno_cols, drop = FALSE]
  rownames(anno) <- anno[[ID]]
  anno[[ID]] <- NULL
  anno[] <- lapply(anno, as.character)

  # Get heatmap palette
  mapal <- palettes(
    category = "heatmap",
    palette = palette_for_heatmape,
    counts = 200,
    show_col = show_col
  )

  # Calculate height
  if (is.null(height)) {
    height <- 2 + length(feas) * 0.25
  }

  # Generate group colors
  cluster_colors <- list()

  # Primary group colors
  mycols1 <- .get_group_colors(cols1, palette1, seed, show_col)
  lev1 <- unique(as.character(input[[group]]))
  lev1 <- lev1[!is.na(lev1)]
  cluster_colors[[group]] <- stats::setNames(mycols1[seq_along(lev1)], lev1)

  # Secondary group colors
  if (!is.null(group2)) {
    if (!group2 %in% colnames(input)) {
      cli::cli_abort("Group2 column {.val {group2}} not found in input")
    }
    mycols2 <- .get_group_colors(cols2, palette2, seed, show_col)
    lev2 <- unique(as.character(input[[group2]]))
    lev2 <- lev2[!is.na(lev2)]
    cluster_colors[[group2]] <- stats::setNames(mycols2[seq_along(lev2)], lev2)
  }

  # Tertiary group colors
  if (!is.null(group3)) {
    if (!group3 %in% colnames(input)) {
      cli::cli_abort("Group3 column {.val {group3}} not found in input")
    }
    mycols3 <- .get_group_colors(cols3, palette3, seed, show_col)
    lev3 <- unique(as.character(input[[group3]]))
    lev3 <- lev3[!is.na(lev3)]
    cluster_colors[[group3]] <- stats::setNames(mycols3[seq_along(lev3)], lev3)
  }

  print(cluster_colors)

  # Create heatmap
  mat <- as.matrix(eset)
  col_fun <- grDevices::colorRampPalette(mapal)

  # Column annotation
  ha_top <- NULL
  if (!is.null(anno) && ncol(anno) > 0) {
    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      df = anno,
      col = cluster_colors,
      annotation_name_gp = grid::gpar(fontsize = 6)
    )
  }

  # Draw heatmap
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "value",
    col = col_fun(256),
    cluster_rows = TRUE,
    cluster_columns = cluster_cols,
    show_row_names = TRUE,
    show_column_names = show_colnames,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6, angle = 45),
    top_annotation = ha_top,
    column_title = NULL,
    row_title = NULL
  )

  # Save output
  outfile <- file.path(
    path_obj$folder_name,
    paste0(file_name_prefix, "-pheatmap-", group, ".", fig.type)
  )

  if (fig.type == "pdf") {
    grDevices::pdf(outfile, width = width, height = height)
  } else if (fig.type %in% c("png", "jpg", "jpeg", "tiff")) {
    # Handle other formats if needed
    cli::cli_alert_warning("Only PDF format is fully supported; saving as PDF")
    grDevices::pdf(outfile, width = width, height = height)
  }

  ComplexHeatmap::draw(ht)
  grDevices::dev.off()

  cli::cli_alert_success("Heatmap saved to: {.path {outfile}}")

  list(
    p_anno = anno,
    p_cols = cluster_colors,
    plot = ht,
    eset = eset
  )
}

#' Get Group Colors
#' @keywords internal
#' @noRd
.get_group_colors <- function(cols, palette, seed, show_col) {
  if (length(cols) == 1 && cols %in% c("random", "normal")) {
    mycols <- palettes(
      category = "random",
      palette = palette,
      show_col = FALSE,
      show_message = FALSE
    )

    if (cols == "random") {
      cli::cli_alert_info("Using random seed: {seed}")
      set.seed(seed)
      mycols <- sample(mycols)
      if (show_col) {
        rlang::check_installed("scales")
        scales::show_col(mycols)
      }
    }
    mycols
  } else {
    if (show_col) {
      rlang::check_installed("scales")
      scales::show_col(cols)
    }
    cols
  }
}
