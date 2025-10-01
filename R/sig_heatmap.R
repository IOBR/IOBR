#' Signature Heatmap with Optional Annotations
#'
#' Generates a heatmap of selected features grouped by a categorical variable, with optional conditional (annotation) bars. Supports palette customization, scaling, size controls, and output saving.
#'
#' @param input Data frame containing ID, grouping variable, and feature columns.
#' @param ID Character. Column name for sample identifier. Default "ID".
#' @param features Character vector. Feature (column) names to include.
#' @param group Character. Grouping variable column name.
#' @param condiction Data frame or NULL. Optional annotation table with variable-condition mapping.
#' @param id_condiction Character. Column name in condiction for feature IDs. Default "vars".
#' @param col_condiction Character. Column name in condiction for condition labels. Default "condiction".
#' @param cols_condiction Character vector. Colors for conditions.
#' @param scale Logical. Scale features by column. Default FALSE.
#' @param palette Integer or character. Palette index/name for heatmap colors. Default 2.
#' @param palette_group Character. Palette group for group colors. Default "jama".
#' @param show_col Logical. Whether to display color vector. Default FALSE.
#' @param show_palettes Logical. Whether to print palette options. Default FALSE.
#' @param cols_group Character vector. Custom colors for groups.
#' @param show_plot Logical. Whether to print the heatmap. Default TRUE.
#' @param width Numeric. Plot width in inches. Default 8.
#' @param height Numeric or NULL. Plot height (auto-calculated if NULL).
#' @param size_col Numeric. Font size for column labels. Default 10.
#' @param size_row Numeric. Font size for row labels. Default 8.
#' @param angle_col Numeric. Rotation angle for column labels. Default 90.
#' @param column_title Character or NULL. Column title.
#' @param row_title Character or NULL. Row title.
#' @param show_heatmap_col_name Logical. Show column names. Default FALSE.
#' @param path Character or NULL. Output directory.
#' @param index Integer or NULL. Index appended to filename.
#'
#' @return Invisibly returns path to saved plot; prints heatmap if show_plot = TRUE.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("tcga_stad_sig", package = "IOBR")
#' data("tcga_stad_pdata", package = "IOBR")
#' input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
#' feas <- grep("MCPcounter", colnames(input), value = TRUE)
#' sig_heatmap(input = input, features = feas, group = "subtype", scale = TRUE)
#'
sig_heatmap <- function(input,
                        ID = "ID",
                        features,
                        group,
                        condiction = NULL,
                        id_condiction = "vars",
                        col_condiction = "condiction",
                        cols_condiction = NULL,
                        scale = FALSE,
                        palette = 2,
                        palette_group = "jama",
                        show_col = F,
                        show_palettes = F,
                        cols_group = NULL,
                        show_plot = T,
                        width = 8,
                        height = NULL,
                        size_col = 10,
                        size_row = 8,
                        angle_col = 90,
                        column_title = NULL,
                        row_title = NULL,
                        show_heatmap_col_name = F,
                        path = NULL,
                        index = NULL) {
  if (!is.null(path)) {
    file_store <- path
  } else {
    file_store <- paste0("1-", group, "-relevant-varbiles-heatmap")
  }

  if (!file.exists(file_store)) dir.create(file_store)
  abspath <- paste(getwd(), "/", file_store, "/", sep = "")


  input <- as.data.frame(input)
  features <- features[features %in% colnames(input)]

  colnames(input)[which(colnames(input) == ID)] <- "idd"
  # input<-column_to_rownames(input,var = "ID")

  input <- input[, c("idd", group, features)]
  colnames(input)[which(colnames(input) == group)] <- "target_group"

  input <- input[!is.na(input[, "target_group"]), ]

  pf_long_group <- tidyr::pivot_longer(input, 3:ncol(input), names_to = "variables", values_to = "value")


  if (!is.null(condiction)) {
    condiction <- condiction[, c(id_condiction, col_condiction)]

    colnames(condiction)[which(colnames(condiction) == id_condiction)] <- "vars"
    colnames(condiction)[which(colnames(condiction) == col_condiction)] <- "condiction"

    pf_long_group <- merge(pf_long_group, condiction, by.x = "variables", by.y = "vars", all.x = TRUE, all.y = FALSE)
    pf_long_group$condiction <- ifelse(is.na(pf_long_group$condiction), "NE", pf_long_group$condiction)
    print(head(pf_long_group))
  }
  ###################################################


  if (is.null(height)) {
    height_heatmap <- length(features) * 0.1 + 3
  } else {
    height_heatmap <- height
  }

  ####################################################
  heatmap_col <- palettes(category = "tidyheatmap", palette = palette, show_col = show_col, show_message = show_palettes)


  if (!is.null(cols_group)) {
    color_box <- cols_group
  } else {
    if (is.null(palette)) {
      color_box <- palettes(category = "random", show_col = show_col, show_message = show_palettes)
    } else {
      color_box <- palettes(category = "box", palette = palette_group, show_col = show_col, show_message = show_palettes)
    }
  }

  if (!is.null(condiction)) {
    target_level1 <- unique(as.character(pf_long_group$condiction))
    target_level1 <- target_level1[!is.na(target_level1)]

    n <- length(target_level1)
    # print(n)
    if (is.null(cols_condiction)) {
      color_box1 <- color_box[1:n]
    } else {
      color_box1 <- cols_condiction
    }
  }

  target_level2 <- unique(as.character(pf_long_group$target_group))
  target_level2 <- target_level2[!is.na(target_level2)]
  n <- length(target_level2)
  # print(n)
  color_box2 <- color_box[1:n]
  # print(color_box)
  ####################################################
  if (scale) {
    scale <- "row"
  } else {
    pf_long_group$value[pf_long_group$value > 3] <- 3
    pf_long_group$value[pf_long_group$value < -3] <- -3
    scale <- "none"
  }

  ####################################################
  if (is.null(condiction)) {
    pp <- pf_long_group %>%
      dplyr::group_by(target_group) %>%
      tidyHeatmap::heatmap(
        .column           = idd,
        .row              = variables,
        .value            = value,
        palette_grouping  = list(c(color_box2)),
        scale             = scale,
        column_title      = column_title,
        row_title         = row_title,
        # annotation      = group2,
        palette_value     = heatmap_col,
        show_column_names = show_heatmap_col_name,
        column_names_gp   = grid::gpar(fontsize = size_col),
        row_names_gp      = grid::gpar(fontsize = size_row),
        column_names_rot  = angle_col
      )
  } else {
    pp <- pf_long_group %>%
      dplyr::group_by(target_group, condiction) %>%
      tidyHeatmap::heatmap(
        .column           = idd,
        .row              = variables,
        .value            = value,
        palette_grouping  = list(c(color_box1), c(color_box2)),
        scale             = scale,
        column_title      = column_title,
        row_title         = row_title,
        # annotation      = group2,
        palette_value     = heatmap_col,
        show_column_names = show_heatmap_col_name,
        column_names_gp   = grid::gpar(fontsize = size_col),
        row_names_gp      = grid::gpar(fontsize = size_row),
        column_names_rot  = angle_col
      )
  }


  if (show_plot) print(pp)

  if (is.null(index)) index <- 1
  pp %>% tidyHeatmap::save_pdf(paste0(abspath, index, "-", group, "-tidyheatmap.pdf"),
    width = width,
    height = height_heatmap
  )

  return(pp)
}
