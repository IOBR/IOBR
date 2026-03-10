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
#' @param cols_heatmap Character vector. Custom colors for heatmap (e.g., c("blue", "white", "red")). Default NULL.
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
                        cols_heatmap = NULL, 
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
  # if (!is.null(path)) {
  #   file_store <- path
  # } else {
  #   file_store <- paste0("1-", group, "-relevant-varbiles-heatmap")
  # }
  # 
  # if (!file.exists(file_store)) dir.create(file_store)
  # abspath <- paste(getwd(), "/", file_store, "/", sep = "")
  # 修改（有条件创建）
  if (!is.null(path)) {
    file_store <- path
    if (!file.exists(file_store)) dir.create(file_store)
    abspath <- paste(getwd(), "/", file_store, "/", sep = "")
  } else {
    abspath <- NULL  # 不创建目录
  }


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
  
  
  # heatmap_col <- palettes(category = "tidyheatmap", palette = palette, show_col = show_col, show_message = show_palettes)
  # 
  # # # 转换为 colorRamp2 函数（tidyHeatmap 1.7.0+ 要求）
  # if (length(heatmap_col) >= 3) {
  #   heatmap_col <- circlize::colorRamp2(c(-2, 0, 2), heatmap_col[1:3])
  # } else {
  #   heatmap_col <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  # }
  
  
  # --- 修改开始：支持自定义热图颜色 ---
  # 用户自定义颜色 - 智能映射
  if (!is.null(cols_heatmap)) {
    # 用户自定义颜色 - 智能映射
    n_colors <- length(cols_heatmap)
    
    if (n_colors >= 5) {
      heatmap_col <- circlize::colorRamp2(c(-2, -1, 0, 1, 2), cols_heatmap[1:5])
      message("Using 5-point color mapping")
    } else if (n_colors >= 3) {
      heatmap_col <- circlize::colorRamp2(c(-2, 0, 2), cols_heatmap[1:3])
    } else if (n_colors == 2) {
      message("Only 2 colors provided, using white as midpoint")
      heatmap_col <- circlize::colorRamp2(c(-2, 0, 2), c(cols_heatmap[1], "white", cols_heatmap[2]))
    } else {
      warning("Invalid cols_heatmap, using default")
      heatmap_col <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    }
    
  } else {
    # 默认使用 palettes 函数
    heatmap_col_raw <- palettes(
      category = "tidyheatmap", 
      palette = palette, 
      show_col = show_col, 
      show_message = show_palettes
    )
    
    # 统一处理各种可能的返回值类型
    if (is.function(heatmap_col_raw)) {
      # 如果是函数（如原版的 colorRamp2 对象），直接使用
      heatmap_col <- heatmap_col_raw
    } else {
      # 转换为颜色向量
      if (is.matrix(heatmap_col_raw) && ncol(heatmap_col_raw) == 3) {
        # 从 RGB 矩阵转换
        color_vector <- rgb(heatmap_col_raw[,1], heatmap_col_raw[,2], heatmap_col_raw[,3], 
                            maxColorValue = 1)
      } else if (is.character(heatmap_col_raw)) {
        # 已经是颜色向量
        color_vector <- heatmap_col_raw
      } else {
        # 回退方案
        color_vector <- c("blue", "white", "red")
      }
      
      # 根据颜色数量创建 colorRamp2
      if (length(color_vector) >= 5) {
        heatmap_col <- circlize::colorRamp2(c(-2, -1, 0, 1, 2), color_vector[1:5])
      } else if (length(color_vector) >= 3) {
        heatmap_col <- circlize::colorRamp2(c(-2, 0, 2), color_vector[1:3])
      } else {
        heatmap_col <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
      }
    }
  }
  # --- 修改结束 ---
  

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
  
  # 在 heatmap 调用前，确保数据是标准 data.frame
  pf_long_group <- as.data.frame(pf_long_group)
  
  
  # if (is.null(condiction)) {
  #   pp <- pf_long_group %>%
  #     dplyr::group_by(target_group) %>%
  #     tidyHeatmap::heatmap(
  #       .column           = idd,
  #       .row              = variables,
  #       .value            = value,
  #       palette_grouping  = list(c(color_box2)),
  #       scale             = scale,
  #       column_title      = column_title,
  #       row_title         = row_title,
  #       # annotation      = group2,
  #       palette_value     = heatmap_col,
  #       show_column_names = show_heatmap_col_name,
  #       column_names_gp   = grid::gpar(fontsize = size_col),
  #       row_names_gp      = grid::gpar(fontsize = size_row),
  #       column_names_rot  = angle_col
  #     )
  # } else {
  #   pp <- pf_long_group %>%
  #     dplyr::group_by(target_group, condiction) %>%
  #     tidyHeatmap::heatmap(
  #       .column           = idd,
  #       .row              = variables,
  #       .value            = value,
  #       palette_grouping  = list(c(color_box1), c(color_box2)),
  #       scale             = scale,
  #       column_title      = column_title,
  #       row_title         = row_title,
  #       # annotation      = group2,
  #       palette_value     = heatmap_col,
  #       show_column_names = show_heatmap_col_name,
  #       column_names_gp   = grid::gpar(fontsize = size_col),
  #       row_names_gp      = grid::gpar(fontsize = size_row),
  #       column_names_rot  = angle_col
  #     )
  # }
  # 
  
  #修改heatmap调用
  if (is.null(condiction)) {
    pp <- pf_long_group %>%
      dplyr::group_by(target_group) %>%
      tidyHeatmap::heatmap(
        .column = idd,
        .row = variables,
        .value = value,
        palette_grouping = list(c(color_box2)),
        scale = scale,
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
      dplyr::group_by(target_group) %>%
      tidyHeatmap::heatmap(
        .column = idd,
        .row = variables,
        .value = value,
        palette_grouping = list(c(color_box1), c(color_box2)),
        scale = scale,
        column_title = column_title,
        row_title = row_title,
        palette_value = heatmap_col,
        show_column_names = show_heatmap_col_name,
        column_names_gp = grid::gpar(fontsize = size_col),  
        row_names_gp = grid::gpar(fontsize = size_row),      
        column_names_rot = angle_col
      )
  }
  
  if (show_plot) print(pp)

  # if (is.null(index)) index <- 1
  # pp %>% tidyHeatmap::save_pdf(paste0(abspath, index, "-", group, "-tidyheatmap.pdf"),
  #   width = width,
  #   height = height_heatmap
  # )
  
  ## 可选保存：只有调用者显式给路径时才写文件
  if (!is.null(path)) {
    if (is.null(index)) index <- 1
    pp %>% tidyHeatmap::save_pdf(
      paste0(abspath, index, "-", group, "-tidyheatmap.pdf"),
      width = width,
      height =height_heatmap
    )
  }

  return(pp)
}
