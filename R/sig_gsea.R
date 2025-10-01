#' Perform Gene Set Enrichment Analysis (GSEA)
#'
#' @description
#' Conducts Gene Set Enrichment Analysis to identify significantly enriched gene sets
#' from differential gene expression data. Supports MSigDB gene sets or custom gene
#' signatures, and generates comprehensive visualizations and statistical results.
#'
#' @param deg Data frame containing differential expression results with gene symbols
#'   and log fold changes.
#' @param genesets List of custom gene sets for enrichment analysis. If \code{NULL},
#'   MSigDB gene sets are used based on \code{org} and \code{category}. Default is
#'   \code{NULL}.
#' @param path Character string specifying the directory path for saving results.
#'   Default is \code{"1-GSEA-result"}.
#' @param gene_symbol Character string specifying the column name in \code{deg}
#'   containing gene symbols. Default is \code{"symbol"}.
#' @param logfc Character string specifying the column name in \code{deg} containing
#'   log fold change values. Default is \code{"log2FoldChange"}.
#' @param org Character string specifying the organism. Options are \code{"hsa"}
#'   (Homo sapiens) or \code{"mus"} (Mus musculus). Default is \code{"hsa"}.
#' @param msigdb Logical indicating whether to use MSigDB gene sets. Default is
#'   \code{TRUE}.
#' @param category Character string specifying the MSigDB category (e.g., \code{"H"}
#'   for Hallmark, \code{"C2"} for curated gene sets). Default is \code{"H"}.
#' @param subcategory Character string specifying the MSigDB subcategory to filter
#'   gene sets. Default is \code{NULL}.
#' @param palette_bar Character string or integer specifying the color palette for
#'   bar plots. Default is \code{"jama"}.
#' @param palette_gsea Integer specifying the color palette for GSEA plots. Default
#'   is 2.
#' @param show_bar Integer specifying the number of top enriched gene sets to display
#'   in the bar plot. Default is 10.
#' @param show_col Logical indicating whether to display color names in the bar plot.
#'   Default is \code{FALSE}.
#' @param show_plot Logical indicating whether to display GSEA enrichment plots.
#'   Default is \code{FALSE}.
#' @param show_gsea Integer specifying the number of top significant gene sets for
#'   which to generate GSEA plots. Default is 8.
#' @param show_path_n Integer specifying the number of pathways to display in GSEA
#'   plots. Default is 20.
#' @param plot_single_sig Logical indicating whether to generate separate plots for
#'   each significant gene set. Default is \code{TRUE}.
#' @param project Character string specifying the project name for output files.
#'   Default is \code{"custom_sig"}.
#' @param minGSSize Integer specifying the minimum gene set size for analysis.
#'   Default is 10.
#' @param maxGSSize Integer specifying the maximum gene set size for analysis.
#'   Default is 500.
#' @param verbose Logical indicating whether to display progress messages. Default
#'   is \code{TRUE}.
#' @param seed Logical indicating whether to set a random seed for reproducibility.
#'   Default is \code{FALSE}.
#' @param fig.type Character string specifying the file format for saving plots
#'   (e.g., \code{"pdf"}, \code{"png"}). Default is \code{"pdf"}.
#' @param print_bar Logical indicating whether to save and print the bar plot.
#'   Default is \code{TRUE}.
#'
#' @return List containing:
#' \itemize{
#'   \item \code{up}: Data frame of up-regulated enriched gene sets
#'   \item \code{down}: Data frame of down-regulated enriched gene sets
#'   \item \code{gsea_result}: Complete GSEA result object
#'   \item \code{gsea_plots}: List of GSEA enrichment plots (if generated)
#' }
#'
#' @author Dongqiang Zeng
#' @export
#' @import DESeq2
#' @examples
#' # Load example data
#' data("eset_stad", package = "IOBR")
#' data("stad_group", package = "IOBR")
#' # Perform differential expression analysis
#' deg <- iobr_deg(eset = eset_stad, pdata = stad_group, group_id = "subtype",
#'                 pdata_id = "ID", array = FALSE, method = "DESeq2",
#'                 contrast = c("EBV", "GS"), path = "STAD")
#' # Run GSEA with custom gene sets
#' res <- sig_gsea(deg = deg, genesets = signature_tme)
sig_gsea <- function(deg,
                     genesets = NULL,
                     path = NULL,
                     gene_symbol = "symbol",
                     logfc = "log2FoldChange",
                     org = "hsa",
                     msigdb = TRUE,
                     category = "H",
                     subcategory = NULL,
                     palette_bar = "jama",
                     palette_gsea = 2,
                     show_bar = 10,
                     show_col = FALSE,
                     show_plot = FALSE,
                     show_gsea = 8,
                     show_path_n = 20,
                     plot_single_sig = FALSE,
                     project = "custom_sig",
                     minGSSize = 10,
                     maxGSSize = 500,
                     verbose = TRUE,
                     seed = FALSE,
                     fig.type = "pdf",
                     print_bar = TRUE) {
  # set path to store enrichment analyses result
  if (is.null(path)) {
    file_store <- paste0("1-GSEA-result")
  } else {
    file_store <- path
  }

  if (!file.exists(file_store)) dir.create(file_store)
  abspath <- paste(getwd(), "/", file_store, "/", sep = "")
  #################################################

  # if(!is.null(input)) (load(paste0(file_source,"/",input)))

  if (gene_symbol != "symbol") {
    colnames(deg)[which(colnames(deg) == "symbol")] <- "symbol_subs"
  }

  colnames(deg)[which(colnames(deg) == gene_symbol)] <- "symbol"
  colnames(deg)[which(colnames(deg) == logfc)] <- "logfc"
  ##################################

  message("`>>>--- Parametar org must be one of: hsa or mus`")
  if (org == "hsa") {
    database <- "org.Hs.eg.db"
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
    entrizid <- clusterProfiler::bitr(deg$symbol,
      fromType = "SYMBOL",
      toType = c("GENENAME", "ENTREZID"),
      OrgDb = database
    )
  } else if (org == "mus") {
    database <- "org.Mm.eg.db"
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
    entrizid <- clusterProfiler::bitr(deg$symbol,
      fromType = "SYMBOL",
      toType = c("GENENAME", "ENTREZID"),
      OrgDb = database
    )
  }

  #################################
  deg <- merge(deg, entrizid, by.x = "symbol", by.y = "SYMBOL", all.x = T, all.y = F)
  ##################################
  gene_id_logfc <- deg %>%
    dplyr::select(ENTREZID, logfc) %>%
    dplyr::distinct(ENTREZID, .keep_all = T) %>%
    filter(!is.na(ENTREZID)) %>%
    filter(!is.na(logfc)) %>%
    arrange(desc(logfc))
  genelist <- gene_id_logfc$logfc
  names(genelist) <- gene_id_logfc$ENTREZID
  genelist <- genelist[order(genelist, decreasing = T)]
  #########################################


  if (is.null(genesets)) {
    if (org == "hsa") {
      species <- "Homo sapiens"
    } else if (org == "mus") {
      species <- "Mus musculus"
    }
    ##################################################
    message(">>>---- Categories that can be choosed... ")

    m_df <- msigdbr::msigdbr(species = species)
    a <- m_df %>%
      dplyr::distinct(gs_cat, gs_subcat) %>%
      dplyr::arrange(gs_cat, gs_subcat)
    print(as.data.frame(a))

    ##################################################
    if (is.null(category)) {
      message(">>>--- Category is NULL, default is Hallmark gene sets...")
      category <- "H"
    }
    term2genes <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
    term2genes <- term2genes %>%
      dplyr::select(gs_name, entrez_gene) %>%
      as.data.frame()
  } else {
    # if(org=="mus"){
    #   genesets<-lapply(genesets, function(x) mus2human_gene)
    # }

    term2genes <- output_sig(genesets, file.name = "sig")
    file.remove("sig.csv")
    term2genes <- reshape2::melt(as.matrix(term2genes))
    term2genes <- term2genes[, -1]
    colnames(term2genes) <- c("gs_name", "symbol")
    category <- project


    if (org == "hsa") {
      database <- "org.Hs.eg.db"
      entrizid <- bitr(term2genes$symbol,
        fromType = "SYMBOL",
        toType = c("GENENAME", "ENTREZID"),
        OrgDb = database
      )
    } else if (org == "mus") {
      database <- "org.Mm.eg.db"
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
      entrizid <- bitr(term2genes$symbol,
        fromType = "SYMBOL",
        toType = c("GENENAME", "ENTREZID"),
        OrgDb = database
      )
    }

    #################################
    term2genes <- merge(term2genes, entrizid, by.x = "symbol", by.y = "SYMBOL", all.x = T, all.y = F)

    term2genes <- term2genes[, c("gs_name", "ENTREZID")]
    colnames(term2genes) <- c("gs_name", "entrez_gene")
  }

  ####################################################
  hall_gsea <- clusterProfiler::GSEA(genelist,
    exponent      = 1,
    # nPerm         = 1000,
    eps           = 0,
    minGSSize     = minGSSize,
    maxGSSize     = maxGSSize,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    TERM2GENE     = term2genes,
    TERM2NAME     = NA,
    verbose       = verbose,
    seed          = FALSE,
    by            = "fgsea"
  )

  if (org == "hsa") {
    hall_gsea <- DOSE::setReadable(hall_gsea, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  } else if (org == "mus") {
    hall_gsea <- DOSE::setReadable(hall_gsea, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
  }

  writexl::write_xlsx(as.data.frame(hall_gsea), paste0(abspath, "1-", category, "_GSEA_significant_results.xlsx"))
  # save(hall_gsea,file = paste(abspath,"2-",category,"_GSEA_result.RData",sep = ""))

  ###########################################
  if (dim(hall_gsea)[1] > 0) {
    if (dim(hall_gsea)[1] < show_gsea) {
      paths <- rownames(hall_gsea@result)
    } else {
      paths <- rownames(hall_gsea[1:show_gsea, ])
    }

    print(paste0(">>>--- Most significant gene sets: "))
    print(paste(1:length(paths), paths, collapse = "; "))
    ###############################################

    # gseacol<- palettes(category = "random", palette = palette_gsea, show_col = FALSE, show_message = F)
    gseacol <- get_cols(palette = palette_gsea, show_col = FALSE)
    ###############################################
    GSEAPLOT <- enrichplot::gseaplot2(
      x = hall_gsea,
      geneSetID = paths,
      subplots = c(1, 2),
      color = gseacol[1:length(paths)],
      pvalue_table = TRUE
    )

    ggplot2::ggsave(
      filename = paste0("2-", category, "_Top_", show_gsea, "_GSEA_plot.", fig.type),
      plot = GSEAPLOT, path = file_store,
      width = 11, height = 7, dpi = 300
    )

    if (show_plot) print(GSEAPLOT)
    ###############################################

    if (plot_single_sig) {
      paths <- rownames(hall_gsea@result)

      paths2 <- paths[1:show_path_n]
      paths2 <- paths2[!is.na(paths2)]

      for (i in 1:length(paths2)) {
        single_path <- paths2[i]
        message(paste0(">>>--- Processing signature: ", single_path))


        ###############################################
        gseaplot <- enrichplot::gseaplot2(
          x = hall_gsea,
          geneSetID = single_path,
          subplots = c(1, 2),
          color = gseacol[i],
          pvalue_table = TRUE
        )
        gseaplot <- gseaplot
        # +design_mytheme(axis_angle = 0, hjust = 0.5)

        if (show_plot) print(gseaplot)
        ggplot2::ggsave(
          filename = paste0(i + 4, "-GSEA_plot-", single_path, ".", fig.type), plot = gseaplot,
          path = file_store, width = 11, height = 7.5, dpi = 300
        )
      }
    }
    ################################################
    down_gogo <- hall_gsea[hall_gsea$pvalue < 0.05 & hall_gsea$enrichmentScore < 0, ]
    if (!dim(down_gogo)[1] == 0) down_gogo$group <- -1
    down_gogo <- down_gogo[1:show_path_n, ]
    down_gogo <- down_gogo[!is.na(down_gogo$p.adjust), ]

    ################################################
    up_gogo <- hall_gsea[hall_gsea$pvalue < 0.05 & hall_gsea$enrichmentScore > 0, ]
    if (!dim(up_gogo)[1] == 0) up_gogo$group <- 1
    up_gogo <- up_gogo[1:show_path_n, ]
    up_gogo <- up_gogo[!is.na(up_gogo$p.adjust), ]


    if (print_bar) {
      ################################################

      gsea_bar <- enrichment_barplot(
        up_terms = up_gogo,
        down_terms = down_gogo,
        palette = palette_bar,
        title = "GSEA-Enrichment",
        width_wrap = 30
      )

      n_bar <- c(dim(up_gogo)[1] + dim(down_gogo)[1])
      height_bar <- 0.5 * n_bar + 3
      ggplot2::ggsave(
        filename = paste0("3-", category, "_GSEA_barplot.", fig.type), plot = gsea_bar,
        path = file_store, width = 6, height = height_bar
      )

      if (show_plot) print(gsea_bar)
      ######################################################
    }

    hall_gsea <- as.tibble(hall_gsea)
    hall_gsea_result <- list(
      up = hall_gsea[hall_gsea$pvalue < 0.05 & hall_gsea$enrichmentScore > 0, ],
      down = hall_gsea[hall_gsea$pvalue < 0.05 & hall_gsea$enrichmentScore < 0, ],
      all = hall_gsea,
      plot_top = GSEAPLOT
    )
    save(hall_gsea_result, file = paste(abspath, "0-", category, "_combined_GSEA_result.RData", sep = ""))
    #######################################################
  } else {
    hall_gsea_result <- hall_gsea
    print(paste0(">>>--- No terms enriched under specific pvalue Cutoff = 0.05"))
  }


  return(hall_gsea_result)
}
