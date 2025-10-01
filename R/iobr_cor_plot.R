#' Integrative Correlation Analysis Between Phenotype and Features
#'
#' @description
#' Performs comprehensive correlation analysis between phenotype data and feature data,
#' supporting both continuous and categorical phenotypes. Filters features based on
#' statistical significance and generates publication-ready visualizations including
#' box plots, heatmaps, and correlation plots.
#'
#' @param pdata_group Data frame containing phenotype data with an identifier column.
#' @param id1 Character string specifying the column name in \code{pdata_group} serving
#'   as the sample identifier.
#' @param feature_data Data frame containing feature data with corresponding identifiers.
#' @param id2 Character string specifying the column name in \code{feature_data} serving
#'   as the sample identifier.
#' @param target Character string specifying the target variable column name for
#'   continuous analysis. Default is \code{NULL}.
#' @param group Character string specifying the grouping variable name for categorical
#'   analysis. Default is \code{"group3"}.
#' @param is_target_continuous Logical indicating whether the target variable is
#'   continuous, which affects grouping strategy. Default is \code{TRUE}.
#' @param padj_cutoff Numeric value specifying the adjusted p-value cutoff for filtering
#'   features. Default is 1.
#' @param index Numeric index used for ordering output file names. Default is 1.
#' @param signature_group Character string specifying the grouping variable for
#'   signatures. Options include \code{"sig_group"} for signature grouping or
#'   \code{"signature_collection"}/\code{"signature_tme"} for gene grouping.
#' @param category Character string specifying the data category: \code{"signature"}
#'   or \code{"gene"}.
#' @param ProjectID Character string specifying the project identifier for file naming.
#' @param feature_limit Integer specifying the maximum number of features to display.
#'   Default is 26.
#' @param character_limit Integer specifying the maximum number of characters for
#'   variable labels. Default is 60.
#' @param palette_box Character string or integer specifying the color palette for
#'   box plots.
#' @param palette_corplot Character string or integer specifying the color palette for
#'   correlation plots.
#' @param palette_heatmap Integer specifying the color palette index for heatmaps.
#' @param show_heatmap_col_name Logical indicating whether to display column names on
#'   heatmaps. Default is \code{TRUE}.
#' @param show_col Logical indicating whether to display color codes for palettes.
#'   Default is \code{FALSE}.
#' @param show_plot Logical indicating whether to display plots. Default is \code{TRUE}.
#' @param path Character string specifying the directory path for saving output files.
#'   Default is \code{NULL}.
#' @param discrete_x Numeric threshold for character length beyond which labels are
#'   discretized. Default is 25.
#' @param show_palettes Logical indicating whether to display color palettes. Default
#'   is \code{FALSE}.
#' @param discrete_width Numeric value specifying the width for label wrapping in plots.
#'   Default is 25.
#' @param fig.type Character string specifying the format for saving figures (\code{"pdf"},
#'   \code{"png"}, etc.). Default is \code{"pdf"}.
#' @param cols_box Character vector of specific colors for box plots. Default is
#'   \code{NULL}.
#'
#' @return Depending on configuration, returns ggplot2 objects (box plots, heatmaps,
#'   correlation plots) and/or a data frame containing statistical analysis results.
#'
#' @author Dongqiang Zeng
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @examples
#' # Create example phenotype and feature data
#' pdata_group <- data.frame(ID = 1:100,
#'                           phenotype = sample(c("Type1", "Type2"), 100, replace = TRUE))
#' feature_data <- data.frame(ID = 1:100,
#'                            Feature1 = rnorm(100),
#'                            Feature2 = rnorm(100))
#' # Perform correlation analysis
#' results <- iobr_cor_plot(
#'   pdata_group = pdata_group, feature_data = feature_data,
#'   id1 = "ID", id2 = "ID", target = "Feature1", is_target_continuous = TRUE,
#'   show_plot = TRUE, path = "path/to/save/results"
#' )
#' print(results)
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
                          signature_group = sig_group,
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
  if (!is.null(path)) {
    file_store <- path
  } else {
    if (!is.null(target)) {
      file_store <- paste0(index, "-1-", ProjectID, "-", target, "-relevant-", category)
    } else {
      file_store <- paste0(index, "-1-", ProjectID, "-", group, "-relevant-", category)
    }
  }

  if (!file.exists(file_store)) dir.create(file_store)
  abspath <- paste(getwd(), "/", file_store, "/", sep = "")

  if (is.null(names(signature_group))) signature_group <- list("signature" = signature_group)
  ##################################################
  pdata_group <- as.data.frame(pdata_group)
  colnames(pdata_group)[which(colnames(pdata_group) == id1)] <- "ID"

  if (!is.null(target) & is_target_continuous) pdata_group[, target] <- as.numeric(pdata_group[, target])

  if (is_target_continuous & !"group2" %in% colnames(pdata_group)) {
    pdata_group$group2 <- ifelse(pdata_group[, target] >= mean(pdata_group[, target]), "High", "Low")
  }

  if (is_target_continuous & !"group3" %in% colnames(pdata_group)) {
    q1 <- quantile(pdata_group[, target], probs = 1 / 3)
    q2 <- quantile(pdata_group[, target], probs = 2 / 3)
    pdata_group$group3 <- ifelse(pdata_group[, target] <= q1, "Low", ifelse(pdata_group[, target] >= q2, "High", "Middle"))
  }

  if (!is.null(target) & is_target_continuous) {
    pdata_group <- pdata_group[, c("ID", target, "group2", "group3")]
  } else {
    pdata_group <- pdata_group[, c("ID", group)]
  }

  # pdata_group<-pdata_group[!is.na(pdata_group[,group]),]
  # pdata_group<-pdata_group[!pdata_group[,group]=="NA",]
  #####################################################

  if (category == "gene") {
    feature_data <- log2eset(feature_data)
    check_eset(feature_data)
    # Gene expression data should be transposed
    feature_data <- tibble::rownames_to_column(as.data.frame(t(feature_data)), var = "ID")
  }

  feature_data <- as.data.frame(feature_data)
  # print(feature_data[1:5,1:5])

  if (category == "signature") {
    if (id2 %in% colnames(feature_data)) {
      colnames(feature_data)[which(colnames(feature_data) == id2)] <- "ID"
    } else {
      stop("id2 is not in the column name of feature_data")
    }
  }

  feature_selected <- feature_manipulation(data = feature_data, feature = setdiff(colnames(feature_data), "ID"))

  # print(feature_selected[1:10])

  feature_data <- feature_data[, colnames(feature_data) %in% c("ID", feature_selected)]

  if (!is.null(target)) {
    if (target %in% colnames(feature_data)) {
      feature_data <- feature_data[, -which(colnames(feature_data) == target)]
    }
  }

  group_list <- signature_group
  panel <- names(signature_group)
  feature_data <- feature_data[, colnames(feature_data) %in% c("ID", unique(unlist(group_list)))]

  if (length(unique(colnames(feature_data))) < length(colnames(feature_data))) stop("There are duplicate variables in the column name of the feature matrix, please modify.")

  # feature_max <- max(feature_data[,colnames(feature_data)%in%feature_selected])
  # if(category == "gene"& feature_max > 50){
  #
  #   rownames(feature_data)<-NULL
  #   feature_data<-column_to_rownames(feature_data,var = "ID")
  #
  #   feature_data<-rownames_to_column(feature_data,var = "ID")
  #   feature_selected<-feature_manipulation(data = feature_data,feature = setdiff(colnames(feature_data),"ID"))
  #   feature_data<-feature_data[,colnames(feature_data)%in%c("ID",feature_selected)]
  #   feature_data<-as.data.frame(feature_data)
  # }

  pf <- merge(pdata_group, feature_data, by = "ID", all = F)
  scale_begin <- length(colnames(pdata_group)) + 1
  pf[, scale_begin:ncol(pf)] <- scale(pf[, scale_begin:ncol(pf)], center = T, scale = T)

  pf_stat <- pf
  ####################################################################
  if (!class(signature_group) == "list") stop(">>> Signature_group must be a list.")
  ####################################################################

  all_sig <- unique(unlist(group_list))

  if (category == "signature") {
    if (length(colnames(pf)[colnames(pf) %in% all_sig]) <= 0) stop(">>> There is no matching signature in signature matrix or signature_group,
                                                              please add new signature groups into sig_group list or generate another one list object")
  }

  if (category == "gene") {
    if (length(colnames(pf)[colnames(pf) %in% all_sig]) <= 0) stop(">>> There is no matching gene in expression set or signature_group,
                                                              please add new signature groups into signature_collections list or generate another one list object")
  }

  if (category == "signature") {
    title.y <- "Signature score"
    title.x <- "Signatures"
  } else {
    title.y <- "Gene expression"
    title.x <- "Signature genes"
  }

  for (x in 1:length(panel)) {
    index_i <- which(names(group_list) == panel[x])[1]
    group_name <- names(group_list)[index_i]

    # print(paste0(">>>  Processing signature: ", group_name))

    features <- group_list[[index_i]]
    features <- features[features %in% colnames(pf)]
    if (length(features) < 2) next

    print(paste0(">>>  Processing signature: ", group_name))
    #####################################
    # if counts of features is larger than feature limit, choose most significant variables
    if (length(features) > feature_limit) {
      if (is_target_continuous == FALSE) {
        eset <- pf[, colnames(pf) %in% c(group, features)]

        if (group == "group3") eset <- eset[!eset$group3 == "Middle", ]

        res <- batch_wilcoxon(data = eset, target = group, feature = setdiff(colnames(eset), group))
        good_features <- high_var_fea(
          result = res,
          target = "sig_names",
          name_padj = "p.adj",
          padj_cutoff = padj_cutoff,
          name_logfc = "statistic",
          logfc_cutoff = 0,
          n = feature_limit / 2
        )
      } else if (is_target_continuous == TRUE) {
        eset <- pf[, colnames(pf) %in% c(target, features)]
        res <- batch_cor(data = eset, target = target, feature = setdiff(colnames(eset), target), method = "spearman")
        good_features <- high_var_fea(
          result = res,
          target = "sig_names",
          name_padj = "p.adj",
          padj_cutoff = padj_cutoff,
          name_logfc = "statistic",
          logfc_cutoff = 0,
          n = feature_limit / 2
        )
      }
      if (length(good_features) <= 2) {
        good_features <- res$sig_names[1:6]
        message(paste0("In panel ", group_name, " : No feature is statistically significant."))
      }

      features <- features[features %in% good_features]
    }

    if (!is.null(target)) {
      pf_inter <- as_tibble(pf[, c(setdiff(colnames(pdata_group), target), target, features)])

      if (target %in% colnames(pf_inter)) index_i_target <- which(colnames(pf_inter) == target)

      pf_long <- tidyr::pivot_longer(pf_inter, c(index_i_target + 1):ncol(pf_inter),
        names_to = "variables", values_to = "value"
      )
    } else {
      pf_inter <- as_tibble(pf[, c("ID", group, features)])
      pf_long <- tidyr::pivot_longer(pf_inter, 3:ncol(pf_inter),
        names_to = "variables", values_to = "value"
      )
    }

    # pf_long[-grep(pf_long$variables,pattern = target),]

    patterns <- tolower(gsub(patterns_to_na, pattern = "\\_", replacement = ""))
    if (tolower(group_name) %in% patterns) {
      pf_long <- remove_names(
        input_df = pf_long,
        variable = "variables",
        patterns_to_na = patterns_to_na,
        patterns_space = c("\\_")
      )
    }


    # pf_long$variables<-gsub(pf_long$variables,pattern = "\\_",replacement = " ")

    pf_long$variables <- substring(pf_long$variables, 1, character_limit)

    if (group == "group3") {
      target_binary <- "group3"
      pf_long_group <- pf_long[!pf_long$group3 == "Middle", ]
    } else if (group == "best_cutoff" & !is.null(target)) {
      target_binary <- paste0(target, "_binary")
      pf_long_group <- pf_long
    } else if (group == "group2") {
      target_binary <- "group2"
      pf_long_group <- pf_long
    } else {
      target_binary <- group
      pf_long_group <- pf_long
    }

    axis_text_size <- 18 - max(nchar(pf_long$variables)) / 7

    pf_long_group_box <- pf_long_group

    if (max(nchar(pf_long_group$variables)) > discrete_x) {
      # pf_long_group_box$nchar_discrete<-ifelse(max(nchar(pf_long_group$variables))> discrete_x,TRUE,FALSE)
      pf_long_group_box$variables <- ifelse(nchar(pf_long_group$variables) > discrete_x, gsub(pf_long_group_box$variables, pattern = "\\_", replacement = " "), pf_long_group_box$variables)
    }


    if (is.null(cols_box)) {
      color_box <- palettes(category = "box", palette = palette_box, show_col = show_col, show_message = show_palettes)
    } else {
      message(">>>-- cols_box must be set if palette_box is NULL...")
      color_box <- cols_box
    }

    ######################################################
    p <- ggboxplot(pf_long_group_box, x = "variables", y = "value", fill = target_binary) +
      scale_fill_manual(values = color_box) +
      ylab(paste0(title.y)) +
      # xlab("")+
      ggtitle(paste0(group_name)) +
      theme_light() +
      theme(
        plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        # axis.text=element_text(size=rel(2.5)),
        axis.text.x = element_text(
          face = "plain", size = axis_text_size,
          angle = 60, hjust = 1, color = "black"
        ), # family="Times New Roman"
        axis.text.y = element_text(face = "plain", size = 15, angle = 0, hjust = 1, color = "black"),
        axis.line = element_line(color = "black", size = .5)
      ) + theme(
        legend.key.size = unit(.3, "inches"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = c(.5, .5),
        legend.box = "horizontal",
        legend.box.just = "top",
        legend.text = element_text(colour = "black", size = 10, face = "plain")
      ) +
      scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = discrete_width))
    #################################################
    max_variables <- max(pf_long_group_box$value)
    group_box <- sym(target_binary)
    pp1 <- p + stat_compare_means(aes(group = !!group_box, label = paste0("p = ", ..p.format..)),
      size = 2.6, label.y = max_variables - 0.3
    )
    pp2 <- p + stat_compare_means(aes(group = !!group_box),
      label = "p.signif",
      size = 6, label.y = max_variables - 0.6
    )
    if (show_plot & length(features) < 13) {
      print(pp1)
    } else if (show_plot & length(features) > 13) {
      print(pp2)
    }
    ###################################################
    plot_width <- length(features) * 0.4 + 3
    plot_height <- 4 + max(nchar(pf_long_group_box$variables)) * 0.05
    ###################################################

    if (fig.type == "pdf") {
      fig.type <- "pdf"
    } else {
      fig.type <- "png"
    }

    if (!is.null(target)) {
      ggsave(pp1,
        filename = paste0(index, "-", x, "-1-", ProjectID, "-", target, "-", group_name, "-pvalue-box.", fig.type),
        width = plot_width, height = plot_height, path = file_store
      )
      ggsave(pp2,
        filename = paste0(index, "-", x, "-2-", ProjectID, "-", target, "-", group_name, "-box.", fig.type),
        width = plot_width, height = plot_height, path = file_store
      )
    } else {
      ggsave(pp1,
        filename = paste0(index, "-", x, "-1-", ProjectID, "-", group, "-", group_name, "-pvalue-box.", fig.type),
        width = plot_width, height = plot_height, path = file_store
      )
      ggsave(pp2,
        filename = paste0(index, "-", x, "-2-", ProjectID, "-", group, "-", group_name, "-box.", fig.type),
        width = plot_width, height = plot_height, path = file_store
      )
    }
    ####################################################

    #---heatmap------------------------

    colnames(pf_long_group)[which(colnames(pf_long_group) == group)] <- "target_group"
    ###################################################
    pf_long_group$value[pf_long_group$value > 2.5] <- 2.5
    pf_long_group$value[pf_long_group$value < -2.5] <- -2.5

    height_heatmap <- length(features) * 0.2 + 3
    ####################################################
    heatmap_col <- palettes(category = "tidyheatmap", palette = palette_heatmap, show_col = show_col, show_message = show_palettes)
    ####################################################
    pp <- pf_long_group %>%
      group_by(target_group) %>%
      tidyHeatmap::heatmap(
        .column = ID,
        .row = variables,
        .value = value,
        palette_grouping = list(c(color_box)),
        # column_title = group_name,
        # annotation = group2,
        palette_value = heatmap_col,
        show_column_names = show_heatmap_col_name
      )

    if (show_plot) print(pp)

    pp %>% tidyHeatmap::save_pdf(paste0(abspath, index, "-", x, "-3-", ProjectID, "-", group, "-", group_name, "-tidyheatmap.pdf"),
      width = 8,
      height = height_heatmap
    )
    ####################################################
    if (is_target_continuous == TRUE & length(group_list[[index_i]]) <= 20) {
      #------corrplot---------------------------------
      pf_cor <- pf[, colnames(pf) %in% c(target, features)]

      if (tolower(group_name) %in% patterns) {
        patterns <- tolower(gsub(patterns_to_na, pattern = "\\_", replacement = ""))
        pf_cor <- remove_names(
          input_df = pf_cor,
          variable = "colnames",
          patterns_to_na = patterns_to_na,
          patterns_space = NULL
        )
      }

      bbcor <- Hmisc::rcorr(as.matrix(pf_cor), type = "spearman")
      bbcor$P[is.na(bbcor$P)] <- 0
      ###################################
      col <- palettes(category = "heatmap3", palette = palette_corplot, show_col = show_col, show_message = show_palettes)

      width_heatmap <- length(group_list[[index_i]]) * 0.75 + 5
      height_heatmap <- length(group_list[[index_i]]) * 0.75 + 4
      #####################################
      pdf(
        file = paste0(
          abspath, index, "-", x, "-4-", ProjectID, "-", group_name,
          "-associated-", category, "-corplot.pdf"
        ),
        width = width_heatmap, height = height_heatmap
      )
      corrplot::corrplot(bbcor$r,
        type = "lower",
        order = "hclust",
        p.mat = bbcor$P,
        sig.level = 0.05,
        tl.srt = 45,
        tl.col = "black",
        tl.cex = 1.3,
        addrect = 2,
        rect.col = "black",
        rect.lwd = 3,
        col = colorRampPalette(col)(50)
      )
      dev.off()
      corrplot::corrplot(bbcor$r,
        type = "lower",
        order = "hclust",
        p.mat = bbcor$P,
        sig.level = 0.05,
        tl.srt = 45,
        tl.col = "black",
        tl.cex = 1,
        addrect = 2,
        rect.col = "black",
        rect.lwd = 3,
        col = colorRampPalette(col)(50)
      )
      ########################################
      lab_size <- 13 - max(nchar(pf_long_group$variables)) / 4 # size of coefficient
      tl_cex <- 20 - max(nchar(pf_long_group$variables)) / 9 # size of signature name
      p <- ggcorrplot::ggcorrplot(bbcor$r,
        hc.order = TRUE, type = "lower", p.mat = bbcor$P, lab = TRUE,
        pch.cex = 4.3,
        lab_size = lab_size,
        tl.cex = tl_cex,
        title = names(group_list)[index_i],
        ggtheme = ggplot2::theme_bw,
        colors = col
      ) +
        theme(plot.title = element_text(size = rel(2.5), hjust = 0.5))
      ######################################
      ggsave(p,
        filename = paste0(
          index, "-", x, "-5-", ProjectID, "-", group_name,
          "-associated-", category, "-corplot.", fig.type
        ),
        width = 12, height = 12.8,
        path = file_store
      )
    }
  }


  if (!is_target_continuous) {
    if (group == "group3") {
      pf_stat <- pf_stat[!pf_stat$group3 == "Middle", ]
    }
    if (length(unique(pf_stat[, group])) == 2) {
      print(">>> Proportion of two groups:")
      print(summary(as.factor(pf_stat[, group])))
      eset <- pf_stat
      feas <- colnames(pf_stat)[scale_begin:ncol(pf_stat)]
      levels <- c(as.character(unique(pf_stat[, group])))
      res <- batch_wilcoxon(data = eset, target = group, feature = feas, feature_manipulation = T)

      res <- tibble::as_tibble(res)
    } else {
      stop("Only two categorical variables can support statistical difference calculation")
    }
  }
  if (is_target_continuous) {
    if (is.null(target)) stop("target must be define...")
    # eset<-pf[,colnames(pf)%in%c(target,setdiff(features,c("ID",target)))]
    eset <- pf_stat
    feas <- colnames(pf_stat)[scale_begin:ncol(pf_stat)]
    res <- batch_cor(data = eset, target = target, feature = feas, method = "spearman")
    res <- tibble::as_tibble(res)
  }
  if (!is.null(res)) {
    return(res)
  }
}
