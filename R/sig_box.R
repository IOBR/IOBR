#' Signature Box Plot with Statistical Comparisons
#'
#' @description
#' Creates box plots to visualize signature distributions across groups with optional
#' statistical pairwise comparisons. Supports both data frames and Seurat objects for
#' single-cell data visualization.
#'
#' @param data Data frame or Seurat object containing the signature and grouping variable.
#' @param signature Character string specifying the column name (or feature name in
#'   Seurat) for the signature values to plot on the y-axis.
#' @param variable Character string specifying the grouping variable column name for
#'   the x-axis.
#' @param palette Character string specifying the color palette name. Default is
#'   \code{"nrc"}.
#' @param cols Character vector of custom fill colors. If \code{NULL}, palette is used.
#'   Default is \code{NULL}.
#' @param jitter Logical indicating whether to add jittered points to the box plot.
#'   Default is \code{FALSE}.
#' @param point_size Numeric value specifying the size of jittered points. Default is 5.
#' @param angle_x_text Numeric value specifying the rotation angle for x-axis labels
#'   (in degrees). Default is 0.
#' @param hjust Numeric value specifying the horizontal justification of x-axis labels.
#'   Default is 0.5.
#' @param show_pvalue Logical indicating whether to display statistical comparison
#'   p-values on the plot. Default is \code{TRUE}.
#' @param return_stat_res Logical indicating whether to return statistical test results
#'   instead of the plot. Default is \code{FALSE}.
#' @param size_of_pvalue Numeric value specifying the font size for p-values. Default
#'   is 6.
#' @param size_of_font Numeric value specifying the base font size. Default is 10.
#' @param assay Character string specifying the assay name (for Seurat objects). Default
#'   is \code{NULL}.
#' @param slot Character string specifying the slot name (for Seurat objects). Default
#'   is \code{"scale.data"}.
#' @param scale Logical indicating whether to scale signature values (z-score
#'   transformation). Default is \code{FALSE}.
#'
#' @return If \code{return_stat_res = FALSE}, returns a ggplot2 object. If
#'   \code{return_stat_res = TRUE}, returns a data frame containing statistical test
#'   results.
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load example data
#' data("tcga_stad_pdata", package = "IOBR")
#' # Create box plot with statistical comparisons
#' sig_box(data = tcga_stad_pdata, signature = "TMEscore_plus",
#'         variable = "subtype", jitter = TRUE, palette = "jco")
sig_box <- function(data, signature, variable, angle_x_text = 0, hjust = 0.5, palette = "nrc", cols = NULL, jitter = FALSE, point_size = 5, size_of_font = 10,
                    size_of_pvalue = 6, show_pvalue = TRUE, return_stat_res = FALSE, assay = NULL, slot = "scale.data", scale = FALSE) {
  if (class(data)[1] == "Seurat") {
    cat(crayon::green(">>>-- Derive matrix data from Seurat object...\n"))
    input <- extract_sc_data(
      sce = data,
      vars = signature,
      assay = assay,
      slot = slot,
      combine_meta_data = TRUE
    )
    data <- input
  }


  data <- as.data.frame(data)
  data <- data[, c(variable, signature)]
  colnames(data)[which(colnames(data) == variable)] <- "variable"
  colnames(data)[which(colnames(data) == signature)] <- "signature"

  if (scale) data[, "signature"] <- as.numeric(scale(data[, "signature"], scale = T, center = T))
  data <- data[!is.na(data$variable), ]
  ####################################
  if (is.null(cols)) {
    cols <- IOBR::palettes(category = "box", palette = palette, show_message = FALSE, show_col = FALSE)
  } else {
    cols <- cols
  }

  p <- ggplot(data, aes(x = variable, y = signature, fill = variable)) +
    geom_boxplot(notch = F, outlier.shape = NULL, outlier.size = 0) +
    scale_fill_manual(values = cols) +
    ylab(signature) +
    xlab(variable)

  comparision <- combn(unique(as.character(data$variable)), 2, simplify = F)

  size_font <- size_of_font * 0.2

  if (angle_x_text %in% c(30, 45, 60)) hjust <- 1

  p <- p + theme_light() +
    theme(
      axis.title.y = element_text(size = rel(size_font)),
      axis.title.x = element_text(size = rel(size_font)),
      axis.text = element_text(size = rel(size_font)),
      axis.text.x = element_text(face = "plain", angle = angle_x_text, hjust = hjust, color = "black"), # family="Times New Roman"

      axis.line = element_line(color = "grey", size = 0.05)
    ) +
    theme(legend.position = "none")


  if (show_pvalue) p <- p + stat_compare_means(comparisons = comparision, size = size_of_pvalue) + stat_compare_means(size = size_of_pvalue)


  res <- compare_means(signature ~ variable, data = data)
  print(res)

  if (jitter) {
    p <- p + geom_point(
      shape = 21,
      size = point_size,
      position = position_jitterdodge(dodge.width = 0.2),
      alpha = .5
    )
  }

  print(p)

  if (return_stat_res) {
    return(res)
  } else {
    return(p)
  }
}
