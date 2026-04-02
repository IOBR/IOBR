#' Batch Signature Box Plots for Group Comparisons
#'
#' @description
#' Generates multiple box plots for specified features (signatures) across groups
#' in the input data. Supports customization of plot appearance, output path,
#' statistical annotation, and compatibility with Seurat objects. Plots are saved
#' to the specified directory or a default folder.
#'
#' @param input Data frame or Seurat object containing the data for analysis.
#' @param vars Character vector. Features or variables to analyze. When
#'   `pattern_vars = TRUE`, these are treated as regular expression patterns.
#' @param groups Character vector. Grouping variable(s) for comparison.
#' @param pattern_vars Logical indicating whether to treat `vars` as regular
#'   expression patterns for matching column names. Default is `FALSE`.
#' @param path Character string or `NULL`. Directory to save plots. If `NULL`,
#'   uses default `"1-sig-box-batch"`.
#' @param index Integer. Starting index for plot filenames. Default is `0`.
#' @param angle_x_text Numeric. Angle of x-axis labels in degrees. Default is `0`.
#' @param hjust Numeric. Horizontal justification of x-axis labels. Default is `0.5`.
#' @param palette Character. Color palette for plots. Default is `"jama"`.
#' @param cols Character vector or `NULL`. Custom colors for plot elements.
#' @param jitter Logical indicating whether to add jittered points. Default is `FALSE`.
#' @param point_size Numeric. Size of points. Default is `5`.
#' @param size_of_font Numeric. Base font size. Default is `8`.
#' @param size_of_pvalue Numeric. Size of p-value text. Default is `4.5`.
#' @param show_pvalue Logical indicating whether to display p-values.
#'   Default is `TRUE`.
#' @param return_stat_res Logical indicating whether to return statistical
#'   results instead of saving plots. Default is `FALSE`.
#' @param assay Character string or `NULL`. Assay type for Seurat objects.
#' @param slot Character. Data slot for Seurat objects. Default is `"scale.data"`.
#' @param scale Logical indicating whether to scale data before analysis.
#'   Default is `FALSE`.
#' @param height Numeric. Height of plots in inches. Default is `5`.
#' @param width Numeric. Width of plots in inches. Default is `3.5`.
#' @param fig_type Character. File format for plots (e.g., `"pdf"`, `"png"`).
#'   Default is `"pdf"`.
#' @param max_count_feas Integer. Maximum number of features to analyze when
#'   `pattern_vars = TRUE`. If matched variables exceed this limit, only the
#'   first `max_count_feas` features are used. Default is `30`.
#'
#' @return If `return_stat_res = TRUE`, returns a data frame of statistical
#'   results; otherwise, invisibly returns the path to saved plots.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' tcga_stad_pdata <- load_data("tcga_stad_pdata")
#' sig_box_batch(
#'   input = tcga_stad_pdata,
#'   vars = c("TMEscore_plus", "GZMB"),
#'   groups = "subtype",
#'   jitter = TRUE,
#'   palette = "jco",
#'   path = tempdir()
#' )
sig_box_batch <- function(input,
                          vars,
                          groups,
                          pattern_vars = FALSE,
                          path = NULL,
                          index = 0,
                          angle_x_text = 0,
                          hjust = 0.5,
                          palette = "jama",
                          cols = NULL,
                          jitter = FALSE,
                          point_size = 5,
                          size_of_font = 8,
                          size_of_pvalue = 4.5,
                          show_pvalue = TRUE,
                          return_stat_res = FALSE,
                          assay = NULL,
                          slot = "scale.data",
                          scale = FALSE,
                          height = 5,
                          width = 3.5,
                          fig_type = "pdf",
                          max_count_feas = 30) {
  # Input validation
  if (!is.data.frame(input) && !inherits(input, "Seurat")) {
    cli::cli_abort("{.arg input} must be a data frame or Seurat object")
  }

  if (!is.character(vars)) {
    cli::cli_abort("{.arg vars} must be a character vector")
  }

  if (!is.character(groups)) {
    cli::cli_abort("{.arg groups} must be a character vector")
  }

  # Handle pattern matching for vars
  if (pattern_vars) {
    rlang::check_installed("stringr")

    vars_matched <- unique(unlist(lapply(vars, function(pattern) {
      colnames(input)[stringr::str_detect(colnames(input), pattern = pattern)]
    })))

    vars_matched <- vars_matched[!is.na(vars_matched)]

    if (length(vars_matched) == 0) {
      cli::cli_abort("No variables matched the specified patterns")
    }

    if (length(vars_matched) > max_count_feas) {
      cli::cli_alert_warning(
        "Matched {length(vars_matched)} variables, using first {max_count_feas}"
      )
      vars <- vars_matched[seq_len(max_count_feas)]
    } else {
      vars <- vars_matched
    }

    cli::cli_alert_info("Variables to analyze: {.val {paste(vars, collapse = ', ')}}")
  }

  # Validate that sig_box_batch is appropriate (needs multiple vars or groups)
  if (length(vars) == 1 && length(groups) == 1) {
    cli::cli_abort(c(
      "{.fn sig_box_batch} requires multiple variables or groups",
      i = "Use {.fn sig_box} for single variable and single group comparisons"
    ))
  }

  # Create output directory
  path <- creat_folder(path %||% "1-sig-box-batch")

  # Process multiple variables
  if (length(vars) > 1) {
    for (i in seq_along(vars)) {
      cli::cli_alert_info("Processing feature: {.val {vars[i]}}")

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
        show_pairwise_p = show_pvalue,
        show_overall_p  = FALSE,
        return_stat_res = return_stat_res,
        assay           = assay,
        slot            = slot,
        scale           = scale
      )

      if (!return_stat_res) {
        ggplot2::ggsave(
          filename = paste0(index + i, "-", vars[i], "-and-", groups, ".", fig_type),
          height   = height,
          width    = width,
          path     = path$folder_name
        )
      }
    }
  }

  # Process multiple groups
  if (length(groups) > 1) {
    base_height <- height

    for (i in seq_along(groups)) {
      cli::cli_alert_info("Processing group: {.val {groups[i]}}")

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
        show_pairwise_p = show_pvalue,
        show_overall_p  = FALSE,
        return_stat_res = return_stat_res,
        assay           = assay,
        slot            = slot,
        scale           = scale
      )

      if (!return_stat_res) {
        # Adjust height based on number of groups
        levs <- unique(input[[groups[i]]])
        levs <- length(levs[!is.na(levs)])
        plot_height <- 2 + base_height * levs

        ggplot2::ggsave(
          filename = paste0(index + i, "-", vars, "-and-", groups[i], ".", fig_type),
          height   = plot_height,
          width    = width,
          path     = path$folder_name
        )
      }
    }
  }

  cli::cli_alert_success("Batch processing complete. Plots saved to: {.path {path$folder_name}}")

  invisible(path$folder_name)
}
