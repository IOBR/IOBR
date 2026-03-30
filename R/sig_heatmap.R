#' Signature Heatmap with Optional Annotations
#'
#' @description
#' Generates a heatmap of selected features grouped by a categorical variable,
#' with optional conditional (annotation) bars. Supports palette customization,
#' scaling, size controls, and output saving.
#'
#' @param input Data frame containing ID, grouping variable, and feature columns.
#' @param id Character string. Column name for sample identifier. Default is `"ID"`.
#' @param features Character vector. Feature (column) names to include in the heatmap.
#' @param group Character string. Grouping variable column name.
#' @param condition Data frame or `NULL`. Optional annotation table with
#'   variable-condition mapping. Default is `NULL`.
#' @param id_condition Character string. Column name in `condition` for feature IDs.
#'   Default is `"vars"`.
#' @param col_condition Character string. Column name in `condition` for condition
#'   labels. Default is `"condition"`.
#' @param cols_condition Character vector. Colors for conditions.
#' @param scale Logical indicating whether to scale values by row. Default is `FALSE`.
#' @param palette Integer or character. Palette index/name for heatmap colors.
#'   Default is `2`.
#' @param cols_heatmap Character vector. Custom colors for heatmap gradient.
#' @param palette_group Character string. Palette name for group colors.
#'   Default is `"jama"`.
#' @param show_col Logical indicating whether to display the color vector.
#'   Default is `FALSE`.
#' @param show_palettes Logical indicating whether to print palette options.
#'   Default is `FALSE`.
#' @param cols_group Character vector. Custom colors for groups.
#' @param show_plot Logical indicating whether to print the heatmap. Default is `TRUE`.
#' @param width Numeric. Plot width in inches. Default is `8`.
#' @param height Numeric or `NULL`. Plot height in inches. Auto-calculated if `NULL`.
#' @param size_col Numeric. Font size for column labels. Default is `10`.
#' @param size_row Numeric. Font size for row labels. Default is `8`.
#' @param angle_col Numeric. Rotation angle for column labels in degrees.
#'   Default is `90`.
#' @param column_title Character string or `NULL`. Title for column annotation.
#' @param row_title Character string or `NULL`. Title for row annotation.
#' @param show_heatmap_col_name Logical indicating whether to show column names.
#'   Default is `FALSE`.
#' @param path Character string or `NULL`. Output directory for saving the heatmap.
#' @param index Integer or `NULL`. Index appended to filename. Default is `NULL`.
#'
#' @return A tidyHeatmap object. Saves PDF only when `path` is provided.
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
#' sig_heatmap(input = input, features = feas, group = "subtype", scale = TRUE)
#' }
sig_heatmap <- function(input,
                        id = "ID",
                        features,
                        group,
                        condition = NULL,
                        id_condition = "vars",
                        col_condition = "condition",
                        cols_condition = NULL,
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
  rlang::check_installed("tidyHeatmap")
  rlang::check_installed("circlize")

  input <- as.data.frame(input)

  # Input validation
  if (!id %in% colnames(input)) {
    cli::cli_abort("ID column {.val {id}} not found in input")
  }

  if (!group %in% colnames(input)) {
    cli::cli_abort("Group column {.val {group}} not found in input")
  }

  features <- unique(features)
  features <- features[features %in% colnames(input)]

  if (length(features) == 0) {
    cli::cli_abort("No valid features found in input")
  }

  cli::cli_alert_info("Creating heatmap with {length(features)} features")

  # Create output directory if needed
  out_dir <- NULL
  if (!is.null(path)) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    out_dir <- normalizePath(path, winslash = "/", mustWork = FALSE)
  }

  # Standardize column names temporarily
  input_copy <- input
  colnames(input_copy)[colnames(input_copy) == id] <- "idd"
  colnames(input_copy)[colnames(input_copy) == group] <- "target_group"

  input_copy <- input_copy[, c("idd", "target_group", features), drop = FALSE]
  input_copy <- input_copy[!is.na(input_copy$target_group), , drop = FALSE]

  # Convert to long format
  pf_long_group <- tidyr::pivot_longer(
    input_copy,
    cols = tidyselect::all_of(features),
    names_to = "variables",
    values_to = "value"
  )

  # Optional row annotation
  if (!is.null(condition)) {
    condition <- as.data.frame(condition)

    if (!all(c(id_condition, col_condition) %in% colnames(condition))) {
      cli::cli_abort(
        "Columns {.val {id_condition}} or {.val {col_condition}} not found in condition"
      )
    }

    condition <- condition[, c(id_condition, col_condition), drop = FALSE]
    colnames(condition) <- c("vars", "condition")

    pf_long_group <- merge(
      pf_long_group,
      condition,
      by.x = "variables",
      by.y = "vars",
      all.x = TRUE,
      all.y = FALSE
    )

    pf_long_group$condition[is.na(pf_long_group$condition)] <- "NE"
  }

  # Calculate plot height
  height_heatmap <- height %||% (length(features) * 0.1 + 3)

  # Build heatmap colors
  heatmap_col <- .build_heatmap_colors(
    cols_heatmap = cols_heatmap,
    palette = palette,
    show_col = show_col,
    show_palettes = show_palettes
  )

  # Get group colors
  color_box <- cols_group %||% palettes(
    category = "box",
    palette = palette_group,
    show_col = show_col,
    show_message = show_palettes
  )

  # Annotation colors for condition
  color_box1 <- NULL
  if (!is.null(condition)) {
    target_level1 <- unique(as.character(pf_long_group$condition))
    target_level1 <- target_level1[!is.na(target_level1)]
    n1 <- length(target_level1)
    color_box1 <- rep(cols_condition %||% color_box, length.out = n1)
  }

  # Group colors for target_group
  target_level2 <- unique(as.character(pf_long_group$target_group))
  target_level2 <- target_level2[!is.na(target_level2)]
  n2 <- length(target_level2)
  color_box2 <- rep(color_box, length.out = n2)

  # Scale or clip values
  if (isTRUE(scale)) {
    scale_mode <- "row"
  } else {
    pf_long_group$value <- pmin(pmax(pf_long_group$value, -3), 3)
    scale_mode <- "none"
  }

  pf_long_group <- as.data.frame(pf_long_group)

  # Build heatmap
  if (is.null(condition)) {
    pp <- pf_long_group %>%
      dplyr::group_by(.data$target_group) %>%
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
      dplyr::group_by(.data$condition, .data$target_group) %>%
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

  # Save to file if path provided
  if (!is.null(out_dir)) {
    idx <- index %||% 1
    outfile <- file.path(out_dir, paste0(idx, "-", group, "-tidyheatmap.pdf"))

    pp %>%
      tidyHeatmap::save_pdf(
        filename = outfile,
        width = width,
        height = height_heatmap
      )

    cli::cli_alert_success("Heatmap saved to: {.path {outfile}}")
  }

  invisible(pp)
}

#' Build Heatmap Colors
#' @keywords internal
#' @noRd
.build_heatmap_colors <- function(cols_heatmap, palette, show_col, show_palettes) {
  if (!is.null(cols_heatmap)) {
    n_colors <- length(cols_heatmap)

    if (n_colors >= 5) {
      return(circlize::colorRamp2(c(-2, -1, 0, 1, 2), cols_heatmap[1:5]))
    } else if (n_colors >= 3) {
      return(circlize::colorRamp2(c(-2, 0, 2), cols_heatmap[1:3]))
    } else if (n_colors == 2) {
      cli::cli_alert_info("Only 2 heatmap colors provided, using white as midpoint")
      return(circlize::colorRamp2(c(-2, 0, 2), c(cols_heatmap[1], "white", cols_heatmap[2])))
    } else {
      cli::cli_warn("Invalid cols_heatmap, using default colors")
      return(circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
    }
  }

  heatmap_col_raw <- palettes(
    category = "tidyheatmap",
    palette = palette,
    show_col = show_col,
    show_message = show_palettes
  )

  if (is.function(heatmap_col_raw)) {
    return(heatmap_col_raw)
  }

  # Convert to color vector
  color_vector <- .convert_to_colors(heatmap_col_raw)

  if (length(color_vector) >= 5) {
    circlize::colorRamp2(c(-2, -1, 0, 1, 2), color_vector[1:5])
  } else if (length(color_vector) >= 3) {
    circlize::colorRamp2(c(-2, 0, 2), color_vector[1:3])
  } else if (length(color_vector) == 2) {
    circlize::colorRamp2(c(-2, 0, 2), c(color_vector[1], "white", color_vector[2]))
  } else {
    circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  }
}

#' Convert various color formats to color vector
#' @keywords internal
#' @noRd
.convert_to_colors <- function(col_raw) {
  if (is.matrix(col_raw) && ncol(col_raw) == 3) {
    max_val <- max(col_raw, na.rm = TRUE)
    grDevices::rgb(
      col_raw[, 1],
      col_raw[, 2],
      col_raw[, 3],
      maxColorValue = ifelse(max_val > 1, 255, 1)
    )
  } else if (is.character(col_raw)) {
    col_raw
  } else {
    c("blue", "white", "red")
  }
}
