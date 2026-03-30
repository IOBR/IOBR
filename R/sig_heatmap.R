#' Signature Heatmap with Optional Annotations
#'
#' Generates a heatmap of selected features grouped by a categorical variable,
#' with optional conditional (annotation) bars. Supports palette customization,
#' scaling, size controls, and output saving.
#'
#' @param input Data frame containing ID, grouping variable, and feature columns.
#' @param ID Character. Column name for sample identifier. Default "ID".
#' @param features Character vector. Feature (column) names to include.
#' @param group Character. Grouping variable column name.
#' @param condiction Data frame or NULL. Optional annotation table with variable-condition mapping.
#' @param id_condiction Character. Column name in condiction for feature IDs. Default "vars".
#' @param col_condiction Character. Column name in condiction for condition labels. Default "condiction".
#' @param cols_condiction Character vector. Colors for conditions.
#' @param scale Logical. Whether to scale values by row. Default FALSE.
#' @param palette Integer or character. Palette index/name for heatmap colors. Default 2.
#' @param cols_heatmap Character vector. Custom colors for heatmap.
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
#' @return A tidyHeatmap object. Saves PDF only when `path` is provided.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' tcga_stad_sig <- load_data("tcga_stad_sig")
#' tcga_stad_pdata <- load_data("tcga_stad_pdata")
#' input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
#' feas <- grep("MCPcounter", colnames(input), value = TRUE)
#' sig_heatmap(input = input, features = feas, group = "subtype", scale = TRUE)
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
                        cols_heatmap = NULL,
                        palette_group = "jama",
                        show_col = FALSE,
                        show_palettes = FALSE,
                        cols_group = NULL,
                        show_plot = TRUE,
                        width = 8,
                        height = NULL,
                        size_col = 10,
                        size_row = 8,
                        angle_col = 90,
                        column_title = NULL,
                        row_title = NULL,
                        show_heatmap_col_name = FALSE,
                        path = NULL,
                        index = NULL) {
  input <- as.data.frame(input)

  if (!ID %in% colnames(input)) {
    stop("`ID` column not found in `input`.")
  }

  if (!group %in% colnames(input)) {
    stop("`group` column not found in `input`.")
  }

  features <- unique(features)
  features <- features[features %in% colnames(input)]

  if (length(features) == 0) {
    stop("No valid `features` found in `input`.")
  }

  ## output path
  if (!is.null(path)) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    out_dir <- normalizePath(path, winslash = "/", mustWork = FALSE)
  } else {
    out_dir <- NULL
  }

  ## standardize names
  colnames(input)[colnames(input) == ID] <- "idd"
  colnames(input)[colnames(input) == group] <- "target_group"

  input <- input[, c("idd", "target_group", features), drop = FALSE]
  input <- input[!is.na(input$target_group), , drop = FALSE]

  ## long format
  pf_long_group <- tidyr::pivot_longer(
    input,
    cols = all_of(features),
    names_to = "variables",
    values_to = "value"
  )

  ## optional row annotation info
  if (!is.null(condiction)) {
    condiction <- as.data.frame(condiction)

    if (!all(c(id_condiction, col_condiction) %in% colnames(condiction))) {
      stop("`id_condiction` or `col_condiction` not found in `condiction`.")
    }

    condiction <- condiction[, c(id_condiction, col_condiction), drop = FALSE]
    colnames(condiction) <- c("vars", "condiction")

    pf_long_group <- merge(
      pf_long_group,
      condiction,
      by.x = "variables",
      by.y = "vars",
      all.x = TRUE,
      all.y = FALSE
    )

    pf_long_group$condiction[is.na(pf_long_group$condiction)] <- "NE"
  }

  ## plot height
  if (is.null(height)) {
    height_heatmap <- length(features) * 0.1 + 3
  } else {
    height_heatmap <- height
  }

  ## heatmap colors
  if (!is.null(cols_heatmap)) {
    n_colors <- length(cols_heatmap)

    if (n_colors >= 5) {
      heatmap_col <- circlize::colorRamp2(
        c(-2, -1, 0, 1, 2),
        cols_heatmap[1:5]
      )
    } else if (n_colors >= 3) {
      heatmap_col <- circlize::colorRamp2(
        c(-2, 0, 2),
        cols_heatmap[1:3]
      )
    } else if (n_colors == 2) {
      message("Only 2 heatmap colors provided, using white as midpoint.")
      heatmap_col <- circlize::colorRamp2(
        c(-2, 0, 2),
        c(cols_heatmap[1], "white", cols_heatmap[2])
      )
    } else {
      warning("Invalid `cols_heatmap`, using default colors.")
      heatmap_col <- circlize::colorRamp2(
        c(-2, 0, 2),
        c("blue", "white", "red")
      )
    }
  } else {
    heatmap_col_raw <- palettes(
      category = "tidyheatmap",
      palette = palette,
      show_col = show_col,
      show_message = show_palettes
    )

    if (is.function(heatmap_col_raw)) {
      heatmap_col <- heatmap_col_raw
    } else {
      if (is.matrix(heatmap_col_raw) && ncol(heatmap_col_raw) == 3) {
        max_val <- max(heatmap_col_raw, na.rm = TRUE)
        color_vector <- grDevices::rgb(
          heatmap_col_raw[, 1],
          heatmap_col_raw[, 2],
          heatmap_col_raw[, 3],
          maxColorValue = ifelse(max_val > 1, 255, 1)
        )
      } else if (is.character(heatmap_col_raw)) {
        color_vector <- heatmap_col_raw
      } else {
        color_vector <- c("blue", "white", "red")
      }

      if (length(color_vector) >= 5) {
        heatmap_col <- circlize::colorRamp2(
          c(-2, -1, 0, 1, 2),
          color_vector[1:5]
        )
      } else if (length(color_vector) >= 3) {
        heatmap_col <- circlize::colorRamp2(
          c(-2, 0, 2),
          color_vector[1:3]
        )
      } else if (length(color_vector) == 2) {
        heatmap_col <- circlize::colorRamp2(
          c(-2, 0, 2),
          c(color_vector[1], "white", color_vector[2])
        )
      } else {
        heatmap_col <- circlize::colorRamp2(
          c(-2, 0, 2),
          c("blue", "white", "red")
        )
      }
    }
  }

  ## grouping colors
  if (!is.null(cols_group)) {
    color_box <- cols_group
  } else {
    color_box <- palettes(
      category = "box",
      palette = palette_group,
      show_col = show_col,
      show_message = show_palettes
    )
  }

  ## annotation colors for condiction
  if (!is.null(condiction)) {
    target_level1 <- unique(as.character(pf_long_group$condiction))
    target_level1 <- target_level1[!is.na(target_level1)]

    n1 <- length(target_level1)

    if (is.null(cols_condiction)) {
      color_box1 <- rep(color_box, length.out = n1)
    } else {
      color_box1 <- rep(cols_condiction, length.out = n1)
    }
  }

  ## grouping colors for target_group
  target_level2 <- unique(as.character(pf_long_group$target_group))
  target_level2 <- target_level2[!is.na(target_level2)]
  n2 <- length(target_level2)
  color_box2 <- rep(color_box, length.out = n2)

  ## scale / clip
  if (isTRUE(scale)) {
    scale_mode <- "row"
  } else {
    pf_long_group$value[pf_long_group$value > 3] <- 3
    pf_long_group$value[pf_long_group$value < -3] <- -3
    scale_mode <- "none"
  }

  pf_long_group <- as.data.frame(pf_long_group)

  ## build heatmap
  if (is.null(condiction)) {
    pp <- pf_long_group %>%
      dplyr::group_by(target_group) %>%
      tidyHeatmap::heatmap(
        .column = idd,
        .row = variables,
        .value = value,
        palette_grouping = list(color_box2),
        scale = scale_mode,
        column_title = column_title,
        row_title = row_title,
        palette_value = heatmap_col,
        show_column_names = show_heatmap_col_name,
        column_names_gp = grid::gpar(fontsize = size_col),
        row_names_gp = grid::gpar(fontsize = size_row),
        column_names_rot = angle_col
      )
  } else {
    pp <- pf_long_group %>%
      dplyr::group_by(condiction, target_group) %>%
      tidyHeatmap::heatmap(
        .column = idd,
        .row = variables,
        .value = value,
        palette_grouping = list(color_box1, color_box2),
        scale = scale_mode,
        column_title = column_title,
        row_title = row_title,
        palette_value = heatmap_col,
        show_column_names = show_heatmap_col_name,
        column_names_gp = grid::gpar(fontsize = size_col),
        row_names_gp = grid::gpar(fontsize = size_row),
        column_names_rot = angle_col
      )
  }

  if (isTRUE(show_plot)) {
    print(pp)
  }

  ## optional save
  if (!is.null(out_dir)) {
    if (is.null(index)) index <- 1
    outfile <- file.path(out_dir, paste0(index, "-", group, "-tidyheatmap.pdf"))

    pp %>%
      tidyHeatmap::save_pdf(
        filename = outfile,
        width = width,
        height = height_heatmap
      )
  }

  return(pp)
}
