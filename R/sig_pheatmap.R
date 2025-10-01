#' Generate Heatmap for Signature Data
#'
#' This function creates a heatmap from signature data with grouping variables,
#' offering flexible options for colors, clustering, and output formats.
#'
#' @param input Data frame with variables in columns.
#' @param feas Vector of feature names (columns) to include in heatmap.
#' @param group Column name for primary grouping variable.
#' @param group2 Optional secondary grouping variable.
#' @param group3 Optional tertiary grouping variable.
#' @param ID Column name for sample identifiers. Default is "ID".
#' @param path Directory to save output files. Default creates "Marker-heatmap-average".
#' @param cols1 Colors for primary group. Default is "random".
#' @param seed Random seed for color generation. Default is 54321.
#' @param show_col Logical indicating whether to display colors. Default is FALSE.
#' @param cluster_cols Logical indicating whether to cluster columns. Default is TRUE.
#' @param palette_for_heatmape Palette number for heatmap. Default is 6.
#' @param scale.matrix Logical indicating whether to scale the matrix. Default is TRUE.
#' @param cellwidth Width of each cell in points. Default is 1.
#' @param cellheight Height of each cell in points. Default is 9.
#' @param fig.type File format for saving. Default is "pdf".
#' @param width Width of saved figure in inches. Default is 6.
#' @param height Height of saved figure in inches. Calculated if NULL.
#' @param file_name_prefix Prefix for saved file name. Default is 1.
#' @param cols2 Colors for secondary group. Default is "random".
#' @param cols3 Colors for tertiary group. Default is "random".
#' @param palette1 Palette for primary group. Default is 1.
#' @param palette2 Palette for secondary group. Default is 2.
#' @param palette3 Palette for tertiary group. Default is 3.
#' @param show_colnames Logical indicating whether to show column names. Default is FALSE.
#'
#' @return A list with annotation data, cluster colors, plot object, and transformed matrix.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("tcga_stad_sig", package = "IOBR")
#' data("tcga_stad_pdata", package = "IOBR")
#' input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
#' feas <- colnames(input)[grep("MCPcounter", colnames(input))]
#' sig_pheatmap(input = input, feas = feas, group = "subtype", scale.matrix = TRUE)
#'
sig_pheatmap <- function(input, feas, group,
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
  if (!is.null(path)) {
    file_store <- path
  } else {
    file_store <- paste0("Marker-heatmap-average")
  }

  path <- creat_folder(file_store)
  ###################################

  input <- as.data.frame(input)
  feas <- feas[feas %in% colnames(input)]

  colnames(input)[which(colnames(input) == ID)] <- "idd"
  ###################################

  eset <- input[, c("idd", feas)]
  rownames(eset) <- NULL
  eset <- column_to_rownames(eset, var = "idd")
  if (scale.matrix == TRUE) eset <- scale(eset, scale = T, center = T)
  eset <- t(eset)
  ##################################

  if (is.null(group3) & is.null(group2)) {
    anno <- input[, c("idd", group)]
  } else if (!is.null(group2) & is.null(goup3)) {
    anno <- input[, c("idd", group, group2)]
  } else {
    anno <- input[, c("idd", group, group2, group3)]
  }

  rownames(anno) <- anno$idd
  anno$idd <- NULL
  ###################################


  mapal <- palettes(category = "heatmap", palette = palette_for_heatmape, counts = 200, show_col = show_col)

  if (is.null(height)) {
    # if(is.null(group)) stop("group must be define")
    height <- 2 + length(feas) * 0.25
  }
  ####################################################


  ##############################################
  if (length(cols1) == 1) {
    if (cols1 == "random") {
      mycols1 <- palettes(category = "random", palette = palette1, show_col = show_col, show_message = FALSE)
      message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
      set.seed(seed)
      mycols1 <- mycols1[sample(length(mycols1), length(mycols1))]
      if (show_col) scales::show_col(mycols1)
    } else if (cols1 == "normal") {
      mycols1 <- palettes(category = "random", palette = palette1, show_col = show_col, show_message = FALSE)
    }
  } else {
    mycols1 <- cols1
    if (show_col) scales::show_col(mycols1)
  }
  ####################################################

  ##############################################
  if (!is.null(group2)) {
    if (length(cols2) == 1) {
      if (cols2 == "random") {
        mycols2 <- palettes(category = "random", palette = palette2, show_col = show_col, show_message = FALSE)
        message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
        set.seed(seed)
        mycols2 <- mycols2[sample(length(mycols2), length(mycols2))]
        if (show_col) scales::show_col(mycols2)
      } else if (cols2 == "normal") {
        mycols2 <- palettes(category = "random", palette = palette2, show_col = show_col, show_message = FALSE)
      }
    } else {
      mycols2 <- cols2
      if (show_col) scales::show_col(mycols2)
    }
  }

  ##############################################
  if (!is.null(group3)) {
    if (length(cols3) == 1) {
      if (cols3 == "random") {
        mycols3 <- palettes(category = "random", palette = palette3, show_col = show_col, show_message = FALSE)
        message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
        set.seed(seed)
        mycols3 <- mycols3[sample(length(mycols3), length(mycols3))]
        if (show_col) scales::show_col(mycols3)
      } else if (cols3 == "normal") {
        mycols3 <- palettes(category = "random", palette = palette3, show_col = show_col, show_message = FALSE)
      }
    } else {
      mycols3 <- cols3
      if (show_col) scales::show_col(mycols3)
    }
  }

  cluster_colors1 <- mycols1[1:length(base::unique(input[, group]))]
  names(cluster_colors1) <- unique(input[, group])
  cluster_colors1 <- list(group = cluster_colors1)
  names(cluster_colors1) <- group
  ######################################################
  if (!is.null(group2)) {
    cluster_colors2 <- mycols2[1:length(base::unique(input[, group2]))]
    names(cluster_colors2) <- unique(input[, group2])
    cluster_colors2 <- list(group2 = cluster_colors2)
    names(cluster_colors2) <- group2
  }
  ######################################################
  if (!is.null(group3)) {
    cluster_colors3 <- mycols3[1:length(base::unique(input[, group3]))]
    names(cluster_colors3) <- unique(input[, group3])
    cluster_colors3 <- list(group3 = cluster_colors3)
    names(cluster_colors3) <- group3
  }
  ######################################################

  if (is.null(group3) & is.null(group2)) {
    cluster_colors <- cluster_colors1
  } else if (!is.null(group2) & is.null(goup3)) {
    cluster_colors <- append(cluster_colors1, cluster_colors2)
  } else {
    cluster_colors <- append(cluster_colors1, cluster_colors2, cluster_colors3)
  }

  print(cluster_colors)
  #######################################################
  # library(pheatmap)
  # pdf(paste0(path$abspath, "2-markers-heatmap-of-average-", group, ".", fig.type), width = width, height = height)
  p <- pheatmap::pheatmap(
    eset,
    color             = mapal, # colorRampPalette(c("darkblue", "white", "red3"))(99), #heatmap color
    scale             = "none",
    cluster_rows      = T,
    cluster_cols      = cluster_cols, # clustered by columns
    cellwidth         = cellwidth,
    cellheight        = cellheight,
    treeheight_col    = 6,
    treeheight_row    = 6,
    clustering_method = "complete",
    show_rownames     = T, # show cluster names
    show_colnames     = show_colnames,
    angle_col         = "45",
    # annotation_row    = annotation_row,
    annotation_col    = anno,
    annotation_colors = cluster_colors,
    fontsize          = 6,
    silent            =  T
  ) # The color for clusters are sames as previous setting
  # dev.off()
  # print(p)
  p
  ggsave(p,
    filename = paste0(file_name_prefix, "-pheatmap-", group, ".", fig.type),
    width = width,
    height = height,
    path = path$folder_name
  )

  res <- list("p_anno" = anno, "p_cols" = cluster_colors, "plot" = p, "eset" = eset)
  return(res)
}
