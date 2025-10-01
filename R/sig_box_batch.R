#' Batch Signature Box Plots for Group Comparisons
#'
#' Generates multiple box plots for specified features (signatures) across groups in the input data. Supports customization of plot appearance, output path, statistical annotation, and compatibility with Seurat objects. Plots are saved to the specified directory or a default folder.
#'
#' @param input Data frame. Input data for analysis.
#' @param vars Character vector. Features or variables to analyze.
#' @param groups Character vector. Grouping variable(s) for comparison.
#' @param path Character or NULL. Directory to save plots. If NULL, uses default "1-sig-box-batch".
#' @param index Integer. Starting index for plot filenames. Default is 0.
#' @param angle_x_text Numeric. Angle of x-axis labels. Default is 0.
#' @param hjust Numeric. Horizontal justification of x-axis labels. Default is 0.5.
#' @param palette Character. Color palette for plots. Default is "jama".
#' @param cols Character vector or NULL. Colors for plot elements.
#' @param jitter Logical. Whether to add jitter to points. Default is FALSE.
#' @param point_size Numeric. Size of points. Default is 5.
#' @param size_of_font Numeric. Font size. Default is 8.
#' @param size_of_pvalue Numeric. Size of p-value text. Default is 4.5.
#' @param show_pvalue Logical. Whether to display p-values. Default is TRUE.
#' @param return_stat_res Logical. Whether to return statistical results instead of saving plots. Default is FALSE.
#' @param assay Character or NULL. Assay type for Seurat objects.
#' @param slot Character. Data slot for Seurat objects. Default is "scale.data".
#' @param scale Logical. Whether to scale data before analysis. Default is FALSE.
#' @param height Numeric. Height of plots. Default is 5.
#' @param width Numeric. Width of plots. Default is 3.5.
#' @param fig_type Character. File format for plots (e.g., "pdf"). Default is "pdf".
#' @param pattern_vars Logical. Whether to treat 'vars' as patterns for matching column names. Default is FALSE.
#'
#' @return If `return_stat_res` is TRUE, returns a data frame of statistical results; otherwise, saves plots to the specified directory.
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data("tcga_stad_pdata", package = "IOBR")
#' sig_box_batch(input = tcga_stad_pdata, vars = c("TMEscore_plus", "GZMB"), groups = "subtype", jitter = TRUE, palette = "jco")
sig_box_batch <- function(input, vars, groups, pattern_vars = FALSE, path = NULL, index = 0, angle_x_text = 0,
                          hjust = 0.5, palette = "jama", cols = NULL, jitter = FALSE, point_size = 5, size_of_font = 8,
                          size_of_pvalue = 4.5, show_pvalue = TRUE, return_stat_res = FALSE, assay = NULL, slot = "scale.data",
                          scale = FALSE, height = 5, width = 3.5, fig_type = "pdf", max_count_feas = 30) {
  if (pattern_vars) {
    vars_com <- c(NULL)
    for (i in 1:length(vars)) {
      vars_loop <- colnames(input)[str_detect(colnames(input), pattern = vars[i])]
      vars_com <- c(vars_com, vars_loop)
    }
    vars <- unique(vars_com[!is.na(vars_com)])
    if (length(vars) > max_count_feas) vars <- vars[1:max_count_feas]
    message(">>>== Variables that will be analysised :")
    message(paste0(vars, collapse = ", "))
  }
  ########################################
  if (length(vars) == 1 & length(groups) == 1) stop(">>>== `sig_box_batch` is suitable for cases where the `vars` or `groups` is greater than one ")

  if (is.null(path)) {
    path <- creat_folder("1-sig-box-batch")
  } else {
    path <- creat_folder(path)
  }
  ########################################

  if (length(vars) > 1) {
    for (i in 1:length(vars)) {
      message(paste0(">>>== Processing feature: ", vars[i], "/n"))
      p <- sig_box(
        data            = input,
        signature       = vars[i],
        variable        = groups,
        angle_x_text    = angle_x_text,
        hjust           = hjust,
        palette         = palette,
        cols            = cols,
        jitter          = jitter,
        point_size      = point_size,
        size_of_font    = size_of_font,
        size_of_pvalue  = size_of_pvalue,
        show_pvalue     = show_pvalue,
        return_stat_res = return_stat_res,
        assay           = assay,
        slot            = slot,
        scale           = scale
      )

      ggsave(filename = paste0(index + i, "-", vars[i], "-and-", groups, ".", fig_type), height = height, width = width, path = path$folder_name)
    }
  }


  if (length(groups) > 1) {
    for (i in 1:length(groups)) {
      message(paste0(">>>== Processing feature: ", groups[i], "/n"))
      p <- sig_box(
        data            = input,
        signature       = vars,
        variable        = groups[i],
        angle_x_text    = angle_x_text,
        hjust           = hjust,
        palette         = palette,
        cols            = cols,
        jitter          = jitter,
        point_size      = point_size,
        size_of_font    = size_of_font,
        size_of_pvalue  = size_of_pvalue,
        show_pvalue     = show_pvalue,
        return_stat_res = return_stat_res,
        assay           = assay,
        slot            = slot,
        scale           = scale
      )

      height_index <- height

      levs <- unique(input[, groups[i]])
      levs <- length(levs[!is.na(levs)])
      height <- 2 + height_index * levs

      ggsave(filename = paste0(index + i, "-", vars, "-and-", groups[i], ".", fig_type), height = height, width = width, path = path$folder_name)
    }
  }

  message(paste0(">>>== Done"))
}
