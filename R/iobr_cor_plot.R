#' Integrative Correlation Analysis Between Phenotype and Features
#'
#' @description
#' Performs comprehensive correlation analysis between phenotype data and feature data,
#' supporting both continuous and categorical phenotypes. Filters features based on
#' statistical significance and generates publication-ready visualizations including
#' box plots, heatmaps, and correlation plots.
#'
#' @param pdata_group Data frame containing phenotype data with an identifier column.
#' @param id1 Character string specifying the column name in `pdata_group` serving
#'   as the sample identifier. Default is `"ID"`.
#' @param feature_data Data frame containing feature data with corresponding identifiers.
#' @param id2 Character string specifying the column name in `feature_data` serving
#'   as the sample identifier. Default is `"ID"`.
#' @param target Character string specifying the target variable column name for
#'   continuous analysis. Default is `NULL`.
#' @param group Character string specifying the grouping variable name for categorical
#'   analysis. Default is `"group3"`.
#' @param is_target_continuous Logical indicating whether the target variable is
#'   continuous, which affects grouping strategy. Default is `TRUE`.
#' @param padj_cutoff Numeric value specifying the adjusted p-value cutoff for filtering
#'   features. Default is `1`.
#' @param index Numeric index used for ordering output file names. Default is `1`.
#' @param signature_group List specifying the grouping variable for signatures.
#'   Options include `"sig_group"` for signature grouping or
#'   `"signature_collection"`/`"signature_tme"` for gene grouping.
#' @param category Character string specifying the data category: `"signature"`
#'   or `"gene"`.
#' @param ProjectID Character string specifying the project identifier for file naming.
#' @param feature_limit Integer specifying the maximum number of features to display.
#'   Default is `26`.
#' @param character_limit Integer specifying the maximum number of characters for
#'   variable labels. Default is `60`.
#' @param palette_box Character string or integer specifying the color palette for
#'   box plots. Default is `"nrc"`.
#' @param palette_corplot Character string or integer specifying the color palette for
#'   correlation plots. Default is `"pheatmap"`.
#' @param palette_heatmap Integer specifying the color palette index for heatmaps.
#'   Default is `2`.
#' @param show_heatmap_col_name Logical indicating whether to display column names on
#'   heatmaps. Default is `FALSE`.
#' @param show_col Logical indicating whether to display color codes for palettes.
#'   Default is `FALSE`.
#' @param show_plot Logical indicating whether to display plots. Default is `FALSE`.
#' @param path Character string specifying the directory path for saving output files.
#'   Default is `NULL`.
#' @param discrete_x Numeric threshold for character length beyond which labels are
#'   discretized. Default is `20`.
#' @param show_palettes Logical indicating whether to display color palettes. Default
#'   is `FALSE`.
#' @param discrete_width Numeric value specifying the width for label wrapping in plots.
#'   Default is `20`.
#' @param fig.type Character string specifying the format for saving figures (`"pdf"`,
#'   `"png"`, etc.). Default is `"pdf"`.
#' @param cols_box Character vector of specific colors for box plots. Default is
#'   `NULL`.
#'
#' @return Depending on configuration, returns ggplot2 objects (box plots, heatmaps,
#'   correlation plots) and/or a data frame containing statistical analysis results.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' pdata_group <- data.frame(
#'   ID = 1:100,
#'   phenotype_score = rnorm(100)
#' )
#'
#' feature_data <- data.frame(
#'   ID = 1:100,
#'   Feature1 = rnorm(100),
#'   Feature2 = rnorm(100),
#'   Feature3 = rnorm(100)
#' )
#'
#' sig_group_example <- list(
#'   signature = c("Feature1", "Feature2", "Feature3")
#' )
#'
#' results <- iobr_cor_plot(
#'   pdata_group = pdata_group,
#'   feature_data = feature_data,
#'   id1 = "ID",
#'   id2 = "ID",
#'   target = "phenotype_score",
#'   is_target_continuous = TRUE,
#'   category = "signature",
#'   signature_group = sig_group_example,
#'   show_plot = FALSE,
#'   path = tempdir()
#' )
#'
#' print(results)
#' }
iobr_cor_plot <- function(pdata_group,
                          id1 = "ID",
                          feature_data,
                          id2 = "ID",
                          target = NULL,
                          group = "group3",
                          is_target_continuous = TRUE,
                          padj_cutoff = 1,
                          index = 1,
                          category = "signature",
                          signature_group = NULL,
                          ProjectID = "TCGA",
                          palette_box = "nrc",
                          cols_box = NULL,
                          palette_corplot = "pheatmap",
                          palette_heatmap = 2,
                          feature_limit = 26,
                          character_limit = 60,
                          show_heatmap_col_name = FALSE,
                          show_col = FALSE,
                          show_plot = FALSE,
                          path = NULL,
                          discrete_x = 20,
                          discrete_width = 20,
                          show_palettes = FALSE,
                          fig.type = "pdf") {
  rlang::check_installed("ggpubr")
  rlang::check_installed("tidyHeatmap")

  if (is.null(signature_group)) {
    cli::cli_abort("{.arg signature_group} must be provided")
  }

  # Create output directory
  file_store <- path %||% if (!is.null(target)) {
    paste0(index, "-1-", ProjectID, "-", target, "-relevant-", category)
  } else {
    paste0(index, "-1-", ProjectID, "-", group, "-relevant-", category)
  }

  if (!dir.exists(file_store)) dir.create(file_store, recursive = TRUE)
  abspath <- file.path(getwd(), file_store, "")

  if (is.null(names(signature_group))) {
    signature_group <- list("signature" = signature_group)
  }

  # Prepare phenotype data
  pdata_group <- as.data.frame(pdata_group)
  colnames(pdata_group)[colnames(pdata_group) == id1] <- "ID"

  if (!is.null(target)) {
    if (!target %in% colnames(pdata_group)) {
      cli::cli_abort("Target {.val {target}} not found in pdata_group")
    }
    if (is_target_continuous) pdata_group[[target]] <- as.numeric(pdata_group[[target]])
  }

  # Create grouping variables
  if (is_target_continuous && !"group2" %in% colnames(pdata_group)) {
    mean_val <- mean(pdata_group[[target]], na.rm = TRUE)
    pdata_group$group2 <- ifelse(pdata_group[[target]] >= mean_val, "High", "Low")
  }

  if (is_target_continuous && !"group3" %in% colnames(pdata_group)) {
    q1 <- stats::quantile(pdata_group[[target]], probs = 1 / 3, na.rm = TRUE)
    q2 <- stats::quantile(pdata_group[[target]], probs = 2 / 3, na.rm = TRUE)
    pdata_group$group3 <- ifelse(
      pdata_group[[target]] <= q1, "Low",
      ifelse(pdata_group[[target]] >= q2, "High", "Middle")
    )
  }

  # Select relevant columns
  if (!is.null(target) && is_target_continuous) {
    pdata_group <- pdata_group[, c("ID", target, "group2", "group3")]
  } else {
    pdata_group <- pdata_group[, c("ID", group)]
  }

  # Prepare feature data
  if (category == "gene") {
    feature_data <- log2eset(feature_data)
    check_eset(feature_data)
    feature_data <- as.data.frame(t(feature_data))
    feature_data <- tibble::rownames_to_column(feature_data, var = "ID")
  }

  feature_data <- as.data.frame(feature_data)

  if (category == "signature") {
    if (!id2 %in% colnames(feature_data)) {
      cli::cli_abort("id2 {.val {id2}} not found in feature_data")
    }
    colnames(feature_data)[colnames(feature_data) == id2] <- "ID"
  }

  feature_selected <- feature_manipulation(
    data = feature_data,
    feature = setdiff(colnames(feature_data), "ID")
  )

  feature_data <- feature_data[, colnames(feature_data) %in% c("ID", feature_selected)]

  if (!is.null(target) && target %in% colnames(feature_data)) {
    feature_data <- feature_data[, colnames(feature_data) != target]
  }

  group_list <- signature_group
  panel <- names(signature_group)
  feature_data <- feature_data[, colnames(feature_data) %in% c("ID", unique(unlist(group_list)))]

  if (any(duplicated(colnames(feature_data)))) {
    cli::cli_abort("Duplicate column names found in feature_data")
  }

  # Merge data
  pf <- merge(pdata_group, feature_data, by = "ID", all = FALSE)
  scale_begin <- length(colnames(pdata_group)) + 1
  pf[, scale_begin:ncol(pf)] <- scale(pf[, scale_begin:ncol(pf)], center = TRUE, scale = TRUE)

  pf_stat <- pf

  # Validate signature_group
  if (!inherits(signature_group, "list")) {
    cli::cli_abort("{.arg signature_group} must be a list")
  }

  all_sig <- unique(unlist(group_list))

  if (length(intersect(colnames(pf), all_sig)) == 0) {
    cli::cli_abort(c(
      "No matching {category} found",
      i = "Check that signature_group matches the data"
    ))
  }

  # Set axis titles
  title.y <- if (category == "signature") "Signature score" else "Gene expression"
  title.x <- if (category == "signature") "Signatures" else "Signature genes"

  # Process each signature group
  for (x in seq_along(panel)) {
    group_name <- panel[x]
    features <- group_list[[x]]
    features <- features[features %in% colnames(pf)]

    if (length(features) < 2) {
      cli::cli_alert_warning("Panel {.val {group_name}} has fewer than 2 features; skipping")
      next
    }

    cli::cli_alert_info("Processing signature: {.val {group_name}}")

    # Select top features if too many
    if (length(features) > feature_limit) {
      if (!is_target_continuous) {
        eset <- pf[, colnames(pf) %in% c(group, features)]
        if (group == "group3") eset <- eset[eset$group3 != "Middle", ]
        res <- batch_wilcoxon(data = eset, target = group, feature = setdiff(colnames(eset), group))
        good_features <- high_var_fea(
          result = res, target = "sig_names", name_padj = "p.adj",
          padj_cutoff = padj_cutoff, name_logfc = "statistic",
          logfc_cutoff = 0, n = feature_limit / 2
        )
      } else {
        eset <- pf[, colnames(pf) %in% c(target, features)]
        res <- batch_cor(
          data = eset, target = target,
          feature = setdiff(colnames(eset), target),
          method = "spearman"
        )
        good_features <- high_var_fea(
          result = res, target = "sig_names", name_padj = "p.adj",
          padj_cutoff = padj_cutoff, name_logfc = "statistic",
          logfc_cutoff = 0, n = feature_limit / 2
        )
      }

      if (length(good_features) <= 2) {
        good_features <- res$sig_names[seq_len(min(6, nrow(res)))]
        cli::cli_alert_warning("Panel {.val {group_name}}: No statistically significant features")
      }

      features <- features[features %in% good_features]
    }

    # Prepare long format data
    if (!is.null(target)) {
      # Select ID, target, group columns, and feature columns
      pf_inter <- tibble::as_tibble(pf[, c("ID", target, group, features)])
      pf_long <- tidyr::pivot_longer(pf_inter, (3 + length(group)):ncol(pf_inter),
        names_to = "variables", values_to = "value"
      )
      pf_long$value <- as.numeric(pf_long$value)
    } else {
      pf_inter <- tibble::as_tibble(pf[, c("ID", group, features)])
      pf_long <- tidyr::pivot_longer(pf_inter, 3:ncol(pf_inter),
        names_to = "variables", values_to = "value",
        values_transform = list(value = as.numeric)
      )
    }

    pf_long$variables <- substring(pf_long$variables, 1, character_limit)

    # Determine target binary variable
    target_binary <- if (group == "group3") {
      pf_long <- pf_long[pf_long$group3 != "Middle", ]
      "group3"
    } else if (group == "group2") {
      "group2"
    } else {
      group
    }

    pf_long_group <- pf_long
    pf_long_group_box <- pf_long_group

    if (max(nchar(as.character(pf_long_group$variables))) > discrete_x) {
      long_idx <- nchar(as.character(pf_long_group_box$variables)) > discrete_x
      pf_long_group_box$variables[long_idx] <- gsub("_", " ", as.character(pf_long_group_box$variables[long_idx]))
    }

    # Get colors
    color_box <- cols_box %||% palettes(
      category = "box", palette = palette_box,
      show_col = show_col, show_message = show_palettes
    )

    # Create box plot
    axis_text_size <- max(8, 18 - max(nchar(as.character(pf_long$variables))) / 7)

    p <- ggpubr::ggboxplot(pf_long_group_box,
      x = "variables", y = "value",
      fill = target_binary
    ) +
      ggplot2::scale_fill_manual(values = color_box) +
      ggplot2::ylab(title.y) +
      ggplot2::ggtitle(group_name) +
      ggplot2::theme_light() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = ggplot2::rel(2), hjust = 0.5),
        axis.title.y = ggplot2::element_text(size = ggplot2::rel(1.5)),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          face = "plain", size = axis_text_size,
          angle = 60, hjust = 1, color = "black"
        ),
        axis.text.y = ggplot2::element_text(
          face = "plain", size = 15, angle = 0,
          hjust = 1, color = "black"
        ),
        axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
        legend.key.size = ggplot2::unit(0.3, "inches"),
        legend.title = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = c(0.5, 0.5),
        legend.box = "horizontal",
        legend.box.just = "top",
        legend.text = ggplot2::element_text(colour = "black", size = 10, face = "plain")
      ) +
      ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = discrete_width))

    max_variables <- max(pf_long_group_box$value, na.rm = TRUE)
    group_box <- rlang::sym(target_binary)

    pp1 <- p + ggpubr::stat_compare_means(
      ggplot2::aes(group = !!group_box, label = paste0("p = ", ggplot2::after_stat(p.format))),
      size = 2.6, label.y = max_variables - 0.3
    )
    pp2 <- p + ggpubr::stat_compare_means(
      ggplot2::aes(group = !!group_box),
      label = "p.signif",
      size = 6, label.y = max_variables - 0.6
    )

    if (show_plot && length(features) < 13) {
      print(pp1)
    } else if (show_plot && length(features) > 13) {
      print(pp2)
    }

    plot_width <- length(features) * 0.4 + 3
    plot_height <- 4 + max(nchar(as.character(pf_long_group_box$variables))) * 0.05

    fig_ext <- if (fig.type == "pdf") "pdf" else "png"

    prefix <- if (!is.null(target)) target else group
    ggplot2::ggsave(pp1,
      filename = paste0(
        index, "-", x, "-1-", ProjectID, "-", prefix,
        "-", group_name, "-pvalue-box.", fig_ext
      ),
      width = plot_width, height = plot_height, path = file_store
    )
    ggplot2::ggsave(pp2,
      filename = paste0(
        index, "-", x, "-2-", ProjectID, "-", prefix,
        "-", group_name, "-box.", fig_ext
      ),
      width = plot_width, height = plot_height, path = file_store
    )

    # Create heatmap
    colnames(pf_long_group)[colnames(pf_long_group) == group] <- "target_group"
    pf_long_group$value <- pmin(pmax(pf_long_group$value, -2.5), 2.5)

    height_heatmap <- length(features) * 0.2 + 3
    heatmap_col <- palettes(
      category = "tidyheatmap", palette = palette_heatmap,
      show_col = show_col, show_message = show_palettes
    )

    pp <- pf_long_group %>%
      dplyr::group_by(.data$target_group) %>%
      tidyHeatmap::heatmap(
        .column = ID, .row = variables, .value = value,
        palette_grouping = list(c(color_box)),
        palette_value = heatmap_col,
        show_column_names = show_heatmap_col_name
      )

    if (show_plot) print(pp)

    pp %>% tidyHeatmap::save_pdf(
      paste0(
        abspath, index, "-", x, "-3-", ProjectID, "-", group, "-", group_name,
        "-tidyheatmap.pdf"
      ),
      width = 8, height = height_heatmap
    )

    # Create correlation plot for continuous targets
    if (is_target_continuous && length(group_list[[x]]) <= 20) {
      pf_cor <- pf[, colnames(pf) %in% c(target, features)]

      rlang::check_installed("Hmisc")
      bbcor <- Hmisc::rcorr(as.matrix(pf_cor), type = "spearman")
      bbcor$P[is.na(bbcor$P)] <- 0

      col <- palettes(
        category = "heatmap3", palette = palette_corplot,
        show_col = show_col, show_message = show_palettes
      )

      width_heatmap <- length(group_list[[x]]) * 0.75 + 5
      height_heatmap <- length(group_list[[x]]) * 0.75 + 4

      grDevices::pdf(
        file = paste0(
          abspath, index, "-", x, "-4-", ProjectID, "-", group_name,
          "-associated-", category, "-corplot.pdf"
        ),
        width = width_heatmap,
        height = height_heatmap
      )
      rlang::check_installed("corrplot")
      corrplot::corrplot(bbcor$r,
        type = "lower", order = "hclust", p.mat = bbcor$P,
        sig.level = 0.05, tl.srt = 45, tl.col = "black",
        tl.cex = 1.3, addrect = 2, rect.col = "black",
        rect.lwd = 3, col = grDevices::colorRampPalette(col)(50)
      )
      grDevices::dev.off()

      corrplot::corrplot(bbcor$r,
        type = "lower", order = "hclust", p.mat = bbcor$P,
        sig.level = 0.05, tl.srt = 45, tl.col = "black",
        tl.cex = 1, addrect = 2, rect.col = "black",
        rect.lwd = 3, col = grDevices::colorRampPalette(col)(50)
      )

      # Create alternative correlation plot with ggplot2
      lab_size <- max(2, 13 - max(nchar(as.character(pf_long_group$variables))) / 4)
      tl_cex <- max(8, 20 - max(nchar(as.character(pf_long_group$variables))) / 9)

      corr <- bbcor$r
      p.mat <- bbcor$P
      corr[upper.tri(corr)] <- NA

      hc <- stats::hclust(stats::dist(corr))
      corr <- corr[hc$order, hc$order]
      p.mat <- p.mat[hc$order, hc$order]

      df <- as.data.frame.table(corr, stringsAsFactors = FALSE)
      colnames(df) <- c("row", "col", "corr")
      df$row <- factor(df$row, levels = rev(unique(df$row)))
      df$col <- factor(df$col, levels = unique(df$col))
      df$stars <- ifelse(p.mat < 0.05 & !is.na(p.mat), "*", "")
      col_fun <- grDevices::colorRampPalette(c("darkblue", "white", "darkred"))

      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$col, y = .data$row,
        fill = .data$corr,
        label = sprintf("%.2f", .data$corr)
      )) +
        ggplot2::geom_tile(color = "grey90") +
        ggplot2::geom_text(color = "black", size = lab_size) +
        ggplot2::geom_text(ggplot2::aes(label = .data$stars),
          color = "red",
          size = lab_size * 1.2, fontface = "bold"
        ) +
        ggplot2::scale_fill_gradientn(colours = col_fun(50), limits = c(-1, 1)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = tl_cex * 0.3),
          axis.text.y = ggplot2::element_text(size = tl_cex * 0.3),
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::ggtitle(group_name)

      ggplot2::ggsave(p,
        filename = paste0(
          index, "-", x, "-5-", ProjectID, "-", group_name,
          "-associated-", category, "-corplot.", fig_ext
        ),
        width = 12, height = 12.8, path = file_store
      )
    }
  }

  # Return statistical results
  if (!is_target_continuous) {
    if (group == "group3") pf_stat <- pf_stat[pf_stat$group3 != "Middle", ]

    if (length(unique(pf_stat[[group]])) == 2) {
      cli::cli_alert_info("Two-group comparison: {table(pf_stat[[group]])}")
      eset <- pf_stat
      feas <- colnames(pf_stat)[scale_begin:ncol(pf_stat)]
      res <- batch_wilcoxon(data = eset, target = group, feature = feas, feature_manipulation = TRUE)
      res <- tibble::as_tibble(res)
    } else {
      cli::cli_abort("Only two categorical variables support statistical difference calculation")
    }
  } else {
    if (is.null(target)) cli::cli_abort("target must be defined for continuous analysis")
    eset <- pf_stat
    feas <- colnames(pf_stat)[scale_begin:ncol(pf_stat)]
    res <- batch_cor(data = eset, target = target, feature = feas, method = "spearman")
    res <- tibble::as_tibble(res)
  }

  invisible(res)
}
