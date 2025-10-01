#' Analyze Mutations Related to Signature Scores
#'
#' This function identifies mutations associated with a specific signature score, performs statistical tests
#' for significance, and generates oncoprints and box plots to visualize relationships.
#'
#' @param mutation_matrix A matrix of mutation data with samples in rows and genes in columns.
#' @param signature_matrix A data frame with sample identifiers and signature scores.
#' @param id_signature_matrix Column name in `signature_matrix` for sample identifiers.
#' @param signature Name of the target signature for analysis.
#' @param min_mut_freq Minimum mutation frequency required for gene inclusion. Default is 0.05.
#' @param plot Logical indicating whether to generate and save plots. Default is TRUE.
#' @param method Statistical test method: "multi" for both Cuzick and Wilcoxon, or "Wilcoxon" only. Default is "multi".
#' @param save_path Directory to save plots and results. If NULL, uses signature name.
#' @param palette Color palette for box plots. Default is "jco".
#' @param show_plot Logical indicating whether to display plots. Default is TRUE.
#' @param show_col Logical indicating whether to show color codes. Default is FALSE.
#' @param width Width of oncoprint plot. Default is 8.
#' @param height Height of oncoprint plot. Default is 4.
#' @param oncoprint_group_by Grouping method for oncoprint: "mean" or "quantile". Default is "mean".
#' @param oncoprint_col Color for mutations in oncoprint. Default is "#224444".
#' @param gene_counts Number of genes to display in oncoprint. Default is 10.
#' @param genes Optional vector of gene names; if NULL, selects based on frequency.
#' @param point_size Size of points in box plot. Default is 4.5.
#' @param point_alpha Transparency of points in box plot. Default is 0.1.
#' @param jitter Logical indicating whether to add jitter to box plot points. Default is FALSE.
#'
#' @return A list containing statistical test results, oncoprint plots, and box plots.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Load mutation and signature data
#' mut_list <- make_mut_matrix(maf = "path_to_maf_file", isTCGA = TRUE, category = "multi")
#' mut <- mut_list$snp
#' results <- find_mutations(mutation_matrix = mut, signature_matrix = tcga_stad_sig,
#'                           id_signature_matrix = "ID", signature = "CD_8_T_effector",
#'                           min_mut_freq = 0.01, plot = TRUE, method = "multi",
#'                           save_path = "path_to_save_results")
find_mutations <- function(mutation_matrix, signature_matrix, id_signature_matrix = "ID", signature,
                           min_mut_freq = 0.05, plot = TRUE, method = "multi", point.alpha = 0.1,
                           save_path = NULL, palette = "jco", show_plot = TRUE,
                           show_col = FALSE, width = 8, height = 4, oncoprint_group_by = "mean",
                           oncoprint_col = "#224444", gene_counts = 10, jitter = FALSE, genes = NULL, point_size = 4.5) {
  if (is.null(save_path)) {
    file_name <- paste0(signature, "-relevant-mutations")
  } else {
    file_name <- save_path
  }

  if (!file.exists(file_name)) dir.create(file_name)
  abspath <- paste0(getwd(), "/", file_name, "/")
  #######################################################
  signature_matrix <- as.data.frame(signature_matrix)

  if (max(mutation_matrix) > 4) {
    mutation_matrix[mutation_matrix >= 3 & mutation_matrix <= 5] <- 3
    mutation_matrix[mutation_matrix > 5] <- 4
  }



  mut2 <- mutation_matrix
  mut2[mut2 >= 1] <- 1
  mut_onco <- mut2

  if (is.null(genes)) {
    mutfreq <- data.frame(head(sort(colSums(mut2), decreasing = T), 500))
    colnames(mutfreq) <- "Freq"
    index <- which(mutfreq$Freq >= dim(mut2)[1] * min_mut_freq)
    index <- max(index)
    input_genes <- names(head(sort(colSums(mut2), decreasing = T), index))
    input_genes <- unique(input_genes)
    input_genes <- input_genes[!is.na(input_genes)]
  } else {
    input_genes <- genes
  }

  mutation_matrix <- mutation_matrix[, colnames(mutation_matrix) %in% input_genes]
  genes <- genes[genes %in% colnames(mutation_matrix)]
  ########################################################

  colnames(signature_matrix)[which(colnames(signature_matrix) == id_signature_matrix)] <- "ID"
  mutation_matrix <- as.data.frame(mutation_matrix)
  mutation_matrix <- tibble::rownames_to_column(mutation_matrix, var = "ID")


  sig_mut <- merge(signature_matrix, mutation_matrix, by = "ID", all = F)
  sig_mut <- tibble::column_to_rownames(sig_mut, var = "ID")


  mytheme <- theme_light() +
    theme(
      plot.title = element_text(size = rel(2.3), hjust = 0.5, face = "italic"),
      axis.title.y = element_text(size = rel(1.8)),
      axis.title.x = element_blank(),
      axis.text.x = element_text(face = "plain", size = 18, angle = 0, color = "black"), # family="Times New Roman"
      axis.text.y = element_text(face = "plain", size = 15, angle = 90, color = "black"), # family="Times New Roman"
      axis.line = element_line(color = "black", size = 0.6)
    ) + theme(
      legend.key.size = unit(.3, "inches"),
      legend.title = element_blank(),
      legend.position = "none",
      legend.direction = "horizontal",
      legend.justification = c(.5, .5),
      legend.box = "vertical",
      legend.box.just = "top",
      legend.text = element_text(colour = "black", size = 10, face = "plain")
    )


  if (method == "multi") {
    if (!is.null(genes)) {
      if (length(genes) < 10) {
        print(genes)
        stop(paste0("Please provide at least 10 genes with mutaion freq larger than ", min_mut_freq))
      }
      input_genes <- genes
    }

    input_genes <- input_genes[input_genes %in% colnames(sig_mut)]

    input <- sig_mut[, c(signature, input_genes)]
    ##############################
    aa <- lapply(input[, input_genes], function(x) PMCMRplus::cuzickTest(input[, 1] ~ x))
    res1 <- data.frame(
      p.value = sapply(aa, getElement, name = "p.value"),
      names = input_genes,
      statistic = sapply(aa, getElement, name = "statistic")
    )
    res1$adjust_pvalue <- p.adjust(res1$p.value, method = "BH", n = length(res1$p.value))
    print(">>>> Result of Cuzick Test")
    res1 <- res1[order(res1$p.value, decreasing = F), ]
    print(res1[1:10, ])

    write.csv(res1, paste0(abspath, "1-cuzickTest-test-relevant-mutations.csv"))


    if (plot) {
      ####################################
      top10_genes <- res1$names[1:10]
      top10_genes <- top10_genes[!is.na(top10_genes)]
      top10_genes <- as.character(top10_genes[top10_genes %in% colnames(input)])
      input_long <- input[, c(signature, top10_genes)]
      input_long <- reshape2::melt(input_long,
        id.vars = 1,
        variable.name = "Gene",
        value.name = "mutation"
      )
      input_long[, signature] <- as.numeric(input_long[, signature])
      input_long$mutation <- as.factor(input_long$mutation)
      ####################################

      pl <- list()
      for (i in 1:length(top10_genes)) {
        gene <- top10_genes[i]
        dd <- input_long[input_long$Gene == gene, ]
        pl[[i]] <- ggplot(dd, aes(x = mutation, y = !!sym(signature), fill = mutation)) +
          geom_boxplot(outlier.shape = NA, outlier.size = -0.5) +
          # geom_jitter(width = 0.25,size= 5.9,alpha=0.75,color ="black")+
          scale_fill_manual(values = palettes(category = "box", palette = palette, show_col = show_col)) +
          mytheme +
          theme(legend.position = "none") +
          ggtitle(paste0(top10_genes[i])) +
          mytheme +
          stat_compare_means(comparisons = combn(as.character(unique(dd[, "mutation"])), 2, simplify = F), size = 6) +
          stat_compare_means(size = 6)

        if (jitter) {
          pl[[i]] <- pl[[i]] + geom_jitter(width = 0.25, size = point_size, alpha = point.alpha, color = "black")
        }
        ggsave(pl[[i]],
          filename = paste0(4 + i, "-1-", gene, "-continue.pdf"),
          width = 4, height = 5.8, path = file_name
        )
      }
      com_plot <- cowplot::plot_grid(pl[[1]], pl[[2]], pl[[3]], pl[[4]], pl[[5]], pl[[6]], pl[[7]], pl[[8]], pl[[9]], pl[[10]],
        labels = "AUTO", ncol = 5, nrow = 2, label_size = 32
      )
      if (show_plot) print(com_plot)
      #####################################
      ggsave(com_plot, filename = "3-Relevant_mutations_Continue.pdf", width = 14, height = 10.5, path = file_name)
      #####################################
    }


    sig_mut2 <- rownames_to_column(sig_mut, var = "ID")
    patr1 <- sig_mut2[, c("ID", signature)]
    part2 <- sig_mut2[, input_genes]
    part2[part2 >= 1] <- 1
    sig_mut2 <- cbind(patr1, part2)
    sig_mut2 <- column_to_rownames(sig_mut2, var = "ID")

    input2 <- sig_mut2

    # print(input2[1:5,1:5])
    aa <- lapply(input2[, input_genes], function(x) wilcox.test(input2[, 1] ~ x))

    res2 <- data.frame(
      p.value = sapply(aa, getElement, name = "p.value"),
      names = input_genes,
      statistic = sapply(aa, getElement, name = "statistic")
    )

    res2$adjust_pvalue <- p.adjust(res2$p.value, method = "BH", n = length(res2$p.value))

    res2 <- res2[order(res2$p.value, decreasing = F), ]
    print(">>> Result of Wilcoxon test (top 10)")
    print(res2[1:10, ])
    write.csv(res2, paste0(abspath, "2-Wilcoxon-test-relevant-mutations.csv"))

    result <- list(
      "cuzick_test" = res1, "wilcoxon_test" = res2,
      "sig_mut_data1" = input, "sig_mut_data2" = input2
    )


    if (plot) {
      ####################################
      top10_genes <- res2$names[1:10]
      top10_genes <- top10_genes[!is.na(top10_genes)]
      top10_genes <- as.character(top10_genes[top10_genes %in% colnames(input2)])

      input_long <- input2[, c(signature, top10_genes)]
      input_long <- reshape2::melt(input_long,
        id.vars = 1,
        variable.name = "Gene",
        value.name = "mutation"
      )
      input_long[, signature] <- as.numeric(input_long[, signature])
      input_long$mutation <- as.factor(input_long$mutation)

      input_long$mutation <- ifelse(input_long$mutation == 0, "WT", "Mutated")
      ####################################

      pl <- list()
      for (i in 1:length(top10_genes)) {
        gene <- top10_genes[i]
        dd <- input_long[input_long$Gene == gene, ]
        pl[[i]] <- ggplot(dd, aes(x = mutation, y = !!sym(signature), fill = mutation)) +
          geom_boxplot(outlier.shape = NA, outlier.size = -0.5) +
          # geom_jitter(width = 0.25,size=5.5,alpha=0.75,color ="black")+
          scale_fill_manual(values = palettes(category = "box", palette = palette, show_col = show_col)) +
          theme(legend.position = "none") +
          ggtitle(paste0(top10_genes[i])) +
          mytheme +
          stat_compare_means(comparisons = combn(as.character(unique(dd[, "mutation"])), 2, simplify = F), size = 6)
        if (jitter) {
          pl[[i]] <- pl[[i]] + geom_jitter(width = 0.25, size = point_size, alpha = point.alpha, color = "black")
        }

        ggsave(pl[[i]],
          filename = paste0(4 + i, "-2-", gene, "-binary.pdf"),
          width = 4, height = 5.8, path = file_name
        )
      }
      com_plot <- cowplot::plot_grid(pl[[1]], pl[[2]], pl[[3]], pl[[4]], pl[[5]], pl[[6]], pl[[7]], pl[[8]], pl[[9]], pl[[10]],
        labels = "AUTO", ncol = 5, nrow = 2, label_size = 32
      )
      if (show_plot) print(com_plot)
      #####################################
      ggsave(com_plot, filename = "4-Relevant_mutations_binary.pdf", width = 14, height = 10.5, path = file_name)
      #####################################
    }
  } else if (tolower(method) == "wilcoxon") {
    if (!is.null(genes)) {
      if (length(genes) < 10) {
        print(genes)
        stop(paste0("Please provide at least 10 genes with mutaion freq larger than ", min_mut_freq))
      }
      input_genes <- genes
    }

    sig_mut2 <- rownames_to_column(sig_mut, var = "ID")

    patr1 <- sig_mut2[, c("ID", signature)]

    input_genes <- input_genes[input_genes %in% colnames(sig_mut2)]

    part2 <- sig_mut2[, input_genes]
    part2[part2 >= 1] <- 1
    sig_mut2 <- cbind(patr1, part2)
    sig_mut2 <- column_to_rownames(sig_mut2, var = "ID")
    input2 <- sig_mut2

    # print(input2[1:5,1:5])
    ##############################
    aa <- lapply(input2[, input_genes], function(x) wilcox.test(input2[, 1] ~ x))

    res <- data.frame(
      p.value = sapply(aa, getElement, name = "p.value"),
      names = input_genes,
      statistic = sapply(aa, getElement, name = "statistic")
    )

    res$adjust_pvalue <- p.adjust(res$p.value, method = "BH", n = length(res$p.value))

    res <- res[order(res$p.value, decreasing = F), ]

    print(">>> Result of Wilcoxon test (top 10): ")
    print(res[1:10, ])

    write.csv(res, paste0(abspath, "0-Wilcoxon-test-relevant-mutations.csv"))

    result <- list("wilcoxon_test" = res, "sig_mut_data" = input2)
    #################################

    if (plot) {
      ####################################
      top10_genes <- res$names[1:10]
      top10_genes <- top10_genes[!is.na(top10_genes)]
      top10_genes <- as.character(top10_genes[top10_genes %in% colnames(input2)])

      input_long <- input2[, c(signature, top10_genes)]
      input_long <- reshape2::melt(input_long,
        id.vars = 1,
        variable.name = "Gene",
        value.name = "mutation"
      )
      input_long[, signature] <- as.numeric(input_long[, signature])
      input_long$mutation <- as.factor(input_long$mutation)

      input_long$mutation <- ifelse(input_long$mutation == 0, "WT", "Mutated")
      ####################################
      pl <- list()
      for (i in 1:length(top10_genes)) {
        gene <- top10_genes[i]
        dd <- input_long[input_long$Gene == gene, ]
        pl[[i]] <- ggplot(dd, aes(x = mutation, y = !!sym(signature), fill = mutation)) +
          geom_boxplot(outlier.shape = NA, outlier.size = NA) +
          # geom_jitter(width = 0.25,size= 3.5,alpha=0.75,color ="black")+
          scale_fill_manual(values = palettes(category = "box", palette = palette, show_col = show_col)) +
          mytheme +
          theme(legend.position = "none") +
          ggtitle(paste0(top10_genes[i])) +
          mytheme +
          stat_compare_means(comparisons = combn(as.character(unique(dd[, "mutation"])), 2, simplify = F), size = 6)

        if (jitter) {
          pl[[i]] <- pl[[i]] + geom_jitter(width = 0.25, size = point_size, alpha = point.alpha, color = "black")
        }

        ggsave(pl[[i]],
          filename = paste0(i, "-1-", gene, "-binary.pdf"),
          width = 4, height = 5.8, path = file_name
        )
      }
      com_plot <- cowplot::plot_grid(pl[[1]], pl[[2]], pl[[3]], pl[[4]], pl[[5]], pl[[6]], pl[[7]], pl[[8]], pl[[9]], pl[[10]],
        labels = "AUTO", ncol = 5, nrow = 2, label_size = 32
      )
      if (show_plot) print(com_plot)
      #####################################
      ggsave(com_plot, filename = "0-Relevant_mutations_binary.pdf", width = 14, height = 10.5, path = file_name)
      #####################################
    }
  }


  # if(!is.null(genes)){
  #   if(length(genes)<10) stop("Please provide at least 10 genes")
  #   genes_for_oncoprint<-genes
  #   genes_for_oncoprint<-genes_for_oncoprint[1:gene_counts]
  #   genes_for_oncoprint<-genes_for_oncoprint[!is.na(genes_for_oncoprint)]
  # }else{
  #   genes_for_oncoprint<-result$wilcoxon_test$names
  #   genes_for_oncoprint<-genes_for_oncoprint[1:gene_counts]
  #   genes_for_oncoprint<-genes_for_oncoprint[!is.na(genes_for_oncoprint)]
  # }
  #
  genes_for_oncoprint <- result$wilcoxon_test$names
  genes_for_oncoprint <- genes_for_oncoprint[1:gene_counts]
  genes_for_oncoprint <- genes_for_oncoprint[!is.na(genes_for_oncoprint)]


  signature_matrix <- signature_matrix[!duplicated(signature_matrix$ID), ]
  signature_matrix <- signature_matrix[signature_matrix$ID %in% rownames(mut_onco), ]

  mut_onco <- mut_onco[rownames(mut_onco) %in% signature_matrix$ID, ]
  ####################################


  if (oncoprint_group_by != "mean" & oncoprint_group_by != "quantile3") {
    pdata_group <- signature_matrix[, c("ID", signature, oncoprint_group_by)]
    pdata_group[, signature] <- as.numeric(pdata_group[, signature])
  } else {
    pdata_group <- signature_matrix[, c("ID", signature)]
    pdata_group[, signature] <- as.numeric(pdata_group[, signature])
  }

  max_sig <- max(pdata_group[, signature], na.rm = TRUE)
  min_sig <- min(pdata_group[, signature], na.rm = TRUE)
  #################################
  if (oncoprint_group_by == "mean") {
    if (!"group2" %in% colnames(pdata_group)) {
      pdata_group$group <- ifelse(pdata_group[, signature] >= mean(pdata_group[, signature]), "High", "Low")
    }
  } else if (oncoprint_group_by == "quantile3") {
    if (!"group3" %in% colnames(pdata_group)) {
      q1 <- quantile(pdata_group[, signature], probs = 1 / 3)
      q2 <- quantile(pdata_group[, signature], probs = 2 / 3)
      pdata_group$group <- ifelse(pdata_group[, signature] <= q1, "Low", ifelse(pdata_group[, signature] >= q2, "High", "Middle"))
    }
  } else {
    if ("group" %in% colnames(pdata_group) & "group" != oncoprint_group_by) stop(" Variable:`group` already exist in colname of pdata!")

    if (!"group" %in% colnames(pdata_group)) {
      colnames(pdata_group)[which(colnames(pdata_group) == oncoprint_group_by)] <- "group"
    }
    # print(head(pdata_group))
    if (nlevels(as.factor(pdata_group$group)) > 2) stop("Levels of `oncoprint_group_by` must less than 3")

    if (!"High" %in% unique(pdata_group$group)) {
      print(summary(as.factor(pdata_group$group)))
      stop("Levels of `oncoprint_group_by` must be `High` or `Low`")
    }
  }


  # print(head(pdata_group))
  idh <- pdata_group[pdata_group$group == "High", "ID"]
  idl <- pdata_group[pdata_group$group == "Low", "ID"]
  pdata1 <- pdata_group[pdata_group$ID %in% idh, ]
  pdata2 <- pdata_group[pdata_group$ID %in% idl, ]

  # library(ComplexHeatmap)
  group_col <- palettes(category = "box", palette = palette, show_col = show_col)

  h1 <- ComplexHeatmap::HeatmapAnnotation(
    Signature_score = anno_barplot(as.numeric(pdata1[, signature]),
      border = FALSE,
      gp = gpar(fill = "#2D004B"),
      axis = TRUE,
      ylim = c(min_sig, max_sig)
    ),
    Group = pdata1$group,
    annotation_height = unit.c(rep(unit(1.5, "cm"), 1), rep(unit(0.5, "cm"), 1)), # unit.c(rep(unit(0.9, "cm"), 5))
    annotation_legend_param = list(
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      ncol = 1
    ),
    gap = unit(c(2, 2), "mm"),
    col = list(Group = c("High" = group_col[1], "Low" = group_col[2])),
    show_annotation_name = TRUE,
    # annotation_name_side="left",
    annotation_name_gp = gpar(fontsize = 12)
  )

  h2 <- ComplexHeatmap::HeatmapAnnotation(
    Signature_score = anno_barplot(as.numeric(pdata2[, signature]),
      border = FALSE,
      gp = gpar(fill = "#2D004B"),
      axis = TRUE,
      ylim = c(min_sig, max_sig)
    ),
    Group = pdata2$group,
    annotation_height = unit.c(rep(unit(1.5, "cm"), 1), rep(unit(0.5, "cm"), 1)), # unit.c(rep(unit(0.9, "cm"), 5))
    annotation_legend_param = list(
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      ncol = 1
    ),
    gap = unit(c(2, 2), "mm"),
    col = list(Group = c("High" = group_col[1], "Low" = group_col[2])),
    show_annotation_name = TRUE,
    # annotation_name_side="left",
    annotation_name_gp = gpar(fontsize = 12)
  )





  col <- c(mut = oncoprint_col)
  mut1 <- t(mut_onco[rownames(mut_onco) %in% idh, colnames(mut_onco) %in% genes_for_oncoprint])
  mut2 <- t(mut_onco[rownames(mut_onco) %in% idl, colnames(mut_onco) %in% genes_for_oncoprint])

  mut1 <- list(mut = mut1)
  mut2 <- list(mut = mut2)

  # width_h<- c(length(idh)/c(length(idh)+length(idl)))*width
  # width_l<- c(length(idl)/c(length(idh)+length(idl)))*width

  #########################################
  ho1 <- ComplexHeatmap::oncoPrint(mut1,
    alter_fun_is_vectorized = FALSE,
    alter_fun = list(mut = function(x, y, w, h) {
      grid.rect(x, y, w * 0.9, h * 0.88,
        gp = gpar(fill = oncoprint_col, col = oncoprint_col)
      )
    }),
    column_title = paste0(" High ", signature),
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(title = "", labels = ""),
    col = col,
    row_names_gp = gpar(fontface = "italic"),
    # width = unit(width_h, "cm"),
    top_annotation = h1
  )


  ho2 <- ComplexHeatmap::oncoPrint(mut2,
    alter_fun_is_vectorized = FALSE,
    alter_fun = list(mut = function(x, y, w, h) {
      grid.rect(x, y, w * 0.9, h * 0.88,
        gp = gpar(fill = oncoprint_col, col = oncoprint_col)
      )
    }),
    column_title = paste0(" Low ", signature),
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(title = "", labels = "Mutation"),
    col = col,
    row_names_gp = gpar(fontface = "italic"),
    # width = unit(width_l, "cm"),
    top_annotation = h2
  )


  p <- ho1 + ho2

  # fig.path<-paste0(getwd(),"/",save_path)
  # save to pdf
  pdf(file.path(abspath, paste0("0-OncoPrint-", signature, ".pdf")), width = width, height = height)
  draw(p)
  invisible(dev.off())
  # print to screen
  # draw(p)
  if (show_plot) {
    print(p)
  }

  result$onco_plot <- p
  result$box_plot <- com_plot
  return(result)
}
