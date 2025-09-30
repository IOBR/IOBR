#' @title sig_heatmap
#' @description The sig_heatmap function is used to generate a heatmap plot based on input data, grouping variables, and optional conditions. The function allows customization of various parameters such as palette selection, scaling, color boxes, plot dimensions, and more. It provides flexibility in visualizing relationships between variables and groups in a concise and informative manner.
#'
#' @param input data frame with `ID`, `group` and features
#' @param ID This parameter specifies the column name in the input data that represents the unique identifier. The default value is "ID".
#' @param features This parameter is a vector of column names representing the features to be used in the heatmap.
#' @param group This parameter is a column name in the input data that represents the grouping variable.
#' @param palette This parameter specifies the palette to be used in the heatmap. It can be a numeric value representing the palette index or a string specifying the palette name. The default value is 2.
#' @param show_col This parameter specifies whether to show the color boxes in the heatmap. If set to TRUE, the color boxes will be shown. The default value is FALSE.
#' @param show_palettes This parameter specifies whether to show the available palettes. If set to TRUE, the available palettes will be shown. The default value is FALSE.
#' @param show_plot This parameter specifies whether to show the resulting heatmap plot. If set to TRUE, the plot will be shown. The default value is TRUE.
#' @param width This parameter specifies the width of the resulting heatmap plot in inches. The default value is 8.
#' @param show_heatmap_col_name This parameter specifies whether to show the column names in the heatmap plot. If set to TRUE, the column names will be shown. The default value is FALSE.
#' @param path This parameter specifies the directory path where the resulting heatmap plot will be saved. If not specified, a default path will be used.
#' @param index default is nullThis parameter specifies the index number to be added to the filename of the resulting heatmap plot. If not specified, the index will be set to 1.
#' @param palette_group This parameter specifies the palette group to be used in the heatmap. It should be a string representing the palette group name. The default value is "jama".
#' @param height This parameter specifies the height of the resulting heatmap plot in inches. If not specified, the height will be calculated based on the number of features.
#' @param size_col This parameter specifies the font size of the column names in the heatmap plot. The default value is 10.
#' @param size_row This parameter specifies the font size of the row names in the heatmap plot. The default value is 8.
#' @param angle_col This parameter specifies the rotation angle of the column names in the heatmap plot. The default value is 90.
#' @param cols_group This parameter is a vector of colors to be used for the grouping variable. If not specified, the colors will be randomly assigned.
#' @param column_title  This parameter specifies the title of the column in the heatmap plot. If not specified, no title will be displayed.
#' @param row_title This parameter specifies the title of the row in the heatmap plot. If not specified, no title will be displayed.
#' @param scale This parameter specifies whether to scale the heatmap by column. If set to TRUE, the heatmap will be scaled by column. The default value is FALSE.
#' @param condiction This parameter is an optional data frame that defines additional conditions to be applied in the heatmap. It should have two columns: one specifying the variables and another specifying the conditions. The default value is NULL.
#' @param id_condiction This parameter specifies the column name in the condiction data frame that represents the variables. The default value is "vars".
#' @param col_condiction This parameter specifies the column name in the condiction data frame that represents the conditions. The default value is "condiction".
#' @param cols_condiction  This parameter is a vector of colors to be used for the condiction variable.
#'
#' @return A heatmap plot object.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("tcga_stad_sig", package = "IOBR")
#' data("tcga_stad_pdata", package = "IOBR")
#' input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
#' feas <- colnames(input)[grep(colnames(input), pattern = "MCPcounter")]
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
