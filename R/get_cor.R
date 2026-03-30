#' Calculate and Visualize Correlation Between Two Variables
#'
#' @description
#' Calculates and visualizes the correlation between two variables with options
#' for scaling, handling missing values, and incorporating grouping data.
#'
#' @param eset Dataset containing the variables (data frame or matrix).
#' @param pdata Optional phenotype data frame. Default is `NULL`.
#' @param var1 Name of the first variable.
#' @param var2 Name of the second variable.
#' @param is.matrix Logical indicating if `eset` is a matrix with features as rows.
#'   Default is `FALSE`.
#' @param id_eset ID column in `eset`. Default is `"ID"`.
#' @param id_pdata ID column in `pdata`. Default is `"ID"`.
#' @param scale Logical indicating whether to scale data. Default is `TRUE`.
#' @param subtype Optional grouping variable for coloring points. Default is `NULL`.
#' @param na.subtype.rm Logical indicating whether to remove NA in subtype.
#'   Default is `FALSE`.
#' @param color_subtype Colors for subtypes. Default is `NULL`.
#' @param palette Color palette name. Default is `"jama"`.
#' @param index Plot index for filename. Default is `NULL` (uses 1).
#' @param method Correlation method: `"spearman"`, `"pearson"`, or `"kendall"`.
#'   Default is `"spearman"`.
#' @param show_cor_result Logical indicating whether to print correlation result.
#'   Default is `TRUE`.
#' @param col_line Color of regression line. Default is `NULL` (auto-determine).
#' @param id Column for point labels. Default is `NULL`.
#' @param show_label Logical indicating whether to show labels. Default is `FALSE`.
#' @param point_size Size of points. Default is 4.
#' @param title Plot title. Default is `NULL`.
#' @param alpha Transparency of points. Default is 0.5.
#' @param title_size Title font size. Default is 1.5.
#' @param text_size Text font size. Default is 10.
#' @param axis_angle Axis label angle. Default is 0.
#' @param hjust Horizontal justification. Default is 0.
#' @param show_plot Logical indicating whether to display plot. Default is `TRUE`.
#' @param save_plot Logical indicating whether to save plot. Default is `FALSE`.
#' @param path Save path. Default is `NULL`.
#' @param fig.format Figure format: `"png"` or `"pdf"`. Default is `"png"`.
#' @param fig.width Figure width in inches. Default is 7.
#' @param fig.height Figure height in inches. Default is 7.3.
#' @param add.hdr.line Logical for adding HDR (high density region) lines.
#'   Default is `FALSE`.
#'
#' @return A ggplot object of the correlation plot.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' eset_tme_stad <- load_data("eset_tme_stad")
#' get_cor(
#'   eset = eset_tme_stad,
#'   is.matrix = TRUE,
#'   var1 = "GZMB",
#'   var2 = "CD274"
#' )
#' }
get_cor <- function(eset,
                    pdata = NULL,
                    var1,
                    var2,
                    is.matrix = FALSE,
                    id_eset = "ID",
                    id_pdata = "ID",
                    scale = TRUE,
                    subtype = NULL,
                    na.subtype.rm = FALSE,
                    color_subtype = NULL,
                    palette = "jama",
                    index = NULL,
                    method = c("spearman", "pearson", "kendall"),
                    show_cor_result = TRUE,
                    col_line = NULL,
                    id = NULL,
                    show_label = FALSE,
                    point_size = 4,
                    title = NULL,
                    alpha = 0.5,
                    title_size = 1.5,
                    text_size = 10,
                    axis_angle = 0,
                    hjust = 0,
                    show_plot = TRUE,
                    save_plot = FALSE,
                    path = NULL,
                    fig.format = "png",
                    fig.width = 7,
                    fig.height = 7.3,
                    add.hdr.line = FALSE) {
  # Validate arguments
  method <- rlang::arg_match(method)
  if (is.null(index)) index <- 1

  # Validate variables
  if (missing(var1) || missing(var2)) {
    cli::cli_abort("Both {.arg var1} and {.arg var2} must be specified.")
  }

  # Prepare data
  data <- .prepare_cor_data(
    eset = eset,
    pdata = pdata,
    var1 = var1,
    var2 = var2,
    is.matrix = is.matrix,
    id_eset = id_eset,
    id_pdata = id_pdata,
    scale = scale
  )

  # Check if variables exist
  if (!var1 %in% colnames(data)) {
    cli::cli_abort("Variable {.val {var1}} not found in data.")
  }
  if (!var2 %in% colnames(data)) {
    cli::cli_abort("Variable {.val {var2}} not found in data.")
  }

  # Remove NA values
  data <- data[!is.na(data[[var1]]) & !is.na(data[[var2]]), , drop = FALSE]

  if (nrow(data) < 3) {
    cli::cli_abort("Insufficient data after removing NA values (need >= 3).")
  }

  cli::cli_alert_info("Calculating {method} correlation (n = {nrow(data)})")

  # Calculate correlation
  cor_result <- stats::cor.test(data[[var1]], data[[var2]], method = method)
  if (show_cor_result) print(cor_result)

  pvalue <- exact_pvalue(data[[var1]], data[[var2]], method = method)
  cli::cli_alert_info("Exact p-value: {format(pvalue, digits = 2, scientific = TRUE)}")

  # Determine line color
  if (is.null(col_line)) {
    col_line <- if (cor_result$estimate > 0) "darkred" else "steelblue"
  }

  # Get color palette
  if (is.null(color_subtype) && !is.null(subtype)) {
    color_subtype <- palettes(
      category = "box",
      palette = palette,
      show_col = FALSE,
      show_message = FALSE
    )
  }

  # Build plot
  p <- .build_cor_plot(
    data = data,
    var1 = var1,
    var2 = var2,
    subtype = subtype,
    na.subtype.rm = na.subtype.rm,
    color_subtype = color_subtype,
    col_line = col_line,
    point_size = point_size,
    alpha = alpha,
    title = title,
    cor_result = cor_result,
    pvalue = pvalue,
    id = id,
    show_label = show_label
  )

  # Apply theme
  theme <- design_mytheme(
    axis_title_size = title_size,
    axis_text_size = text_size,
    hjust = hjust,
    axis_angle = axis_angle
  )

  p <- p + theme +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(
        size = 15,
        hjust = 0,
        face = "italic",
        color = "black"
      )
    )

  # Add HDR lines if requested
  if (add.hdr.line) {
    rlang::check_installed("ggdensity")
    p <- p + ggdensity::geom_hdr_lines()
  }

  if (show_plot) print(p)

  # Save plot if requested
  if (save_plot) {
    .save_cor_plot(
      p = p,
      path = path,
      var1 = var1,
      var2 = var2,
      index = index,
      fig.format = fig.format,
      fig.width = fig.width,
      fig.height = fig.height,
      data = data
    )
  }

  invisible(p)
}

# Helper: Prepare data for correlation
.prepare_cor_data <- function(eset, pdata, var1, var2, is.matrix, id_eset,
                              id_pdata, scale) {
  if (is.null(pdata)) {
    if (is.matrix) {
      if (scale) {
        eset <- scale_matrix(
          matrix = eset,
          log2matrix = FALSE,
          manipulate = FALSE
        )
      }
      data <- as.data.frame(t(eset))
      data$ID <- rownames(data)
    } else {
      data <- as.data.frame(eset)
      if (id_eset %in% colnames(data)) {
        colnames(data)[colnames(data) == id_eset] <- "ID"
      }
      if (scale) {
        numeric_cols <- sapply(data, is.numeric)
        data[, numeric_cols] <- scale(data[, numeric_cols])
      }
    }
  } else {
    # Merge pdata with eset
    pdata <- as.data.frame(pdata)
    colnames(pdata)[colnames(pdata) == id_pdata] <- "ID"

    if (is.matrix) {
      feas <- c(var1, var2)
      feas <- feas[feas %in% rownames(eset)]
      if (length(feas) == 0) {
        cli::cli_abort("No specified variables found in eset.")
      }
      data <- combine_pd_eset(
        eset = eset,
        pdata = pdata,
        id_pdata = "ID",
        feas = feas,
        scale = scale
      )
    } else {
      eset <- as.data.frame(eset)
      if (id_eset %in% colnames(eset)) {
        colnames(eset)[colnames(eset) == id_eset] <- "ID"
      }
      data <- merge(pdata, eset, by = "ID", all = FALSE)
    }
  }

  data
}

# Helper: Build correlation plot
.build_cor_plot <- function(data, var1, var2, subtype, na.subtype.rm,
                            color_subtype, col_line, point_size, alpha,
                            title, cor_result, pvalue, id, show_label) {
  # Base aesthetics
  if (!is.null(subtype)) {
    if (!subtype %in% colnames(data)) {
      cli::cli_warn("Subtype {.val {subtype}} not found. Ignoring.")
      subtype <- NULL
    }
  }

  if (!is.null(subtype)) {
    colnames(data)[colnames(data) == subtype] <- "categorys"

    if (na.subtype.rm) {
      data <- data[!is.na(data$categorys), , drop = FALSE]
    }

    data$categorys <- as.character(data$categorys)
    data$categorys[is.na(data$categorys)] <- "Not_available"
    data$categorys <- as.factor(data$categorys)

    cli::cli_alert_info("Groups: {.val {levels(data$categorys)}}")

    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[var1]],
        y = .data[[var2]],
        colour = .data$categorys
      )
    ) +
      ggplot2::scale_color_manual(values = color_subtype) +
      ggplot2::geom_point(size = point_size, alpha = alpha)
  } else {
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[var1]],
        y = .data[[var2]]
      )
    ) +
      ggplot2::geom_point(
        size = point_size,
        alpha = alpha,
        colour = "black"
      )
  }

  # Add regression line and labels
  p <- p +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = col_line) +
    ggplot2::labs(
      x = var1,
      y = var2,
      title = title,
      subtitle = paste0(
        "r = ", round(as.numeric(cor_result$estimate), 3),
        ", P = ", format(pvalue, digits = 1, scientific = TRUE)
      )
    )

  # Add labels if requested
  if (show_label) {
    if (is.null(id) || !id %in% colnames(data)) {
      cli::cli_warn("Label column not found. Set {.arg id} to show labels.")
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data[[id]]),
        nudge_x = 0.25,
        nudge_y = 0.25,
        check_overlap = TRUE,
        size = 3
      )
    }
  }

  p
}

# Helper: Save correlation plot
.save_cor_plot <- function(p, path, var1, var2, index, fig.format,
                           fig.width, fig.height, data) {
  if (is.null(path)) {
    path <- paste0("1-Cor-of-", var1, "-and-", var2)
  }

  ff <- creat_folder(path)

  filename <- paste0(index, "-", var2, "-", var1, "-correlation.", fig.format)

  ggplot2::ggsave(
    p,
    filename = filename,
    width = fig.width,
    height = fig.height,
    path = ff$folder_name
  )

  save(data, file = paste0(ff$abspath, "0-input-data-", var1, "-", var2, ".RData"))

  cli::cli_alert_success("Plot saved to {.path {ff$folder_name}}")
}
