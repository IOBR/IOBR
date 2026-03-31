#' Perform Gene Set Enrichment Analysis (GSEA)
#'
#' @description
#' Conducts Gene Set Enrichment Analysis to identify significantly enriched gene sets
#' from differential gene expression data. Supports MSigDB gene sets or custom gene
#' signatures, and generates comprehensive visualizations and statistical results.
#'
#' @param deg Data frame containing differential expression results with gene symbols
#'   and log fold changes.
#' @param genesets List of custom gene sets for enrichment analysis. If `NULL`,
#'   MSigDB gene sets are used based on `org` and `category`. Default is `NULL`.
#' @param path Character string specifying the directory path for saving results.
#'   Default is `NULL`.
#' @param gene_symbol Character string specifying the column name in `deg`
#'   containing gene symbols. Default is `"symbol"`.
#' @param logfc Character string specifying the column name in `deg` containing
#'   log fold change values. Default is `"log2FoldChange"`.
#' @param org Character string specifying the organism. Options are `"hsa"`
#'   (Homo sapiens) or `"mus"` (Mus musculus). Default is `"hsa"`.
#' @param msigdb Logical indicating whether to use MSigDB gene sets. Default is
#'   `TRUE`.
#' @param category Character string specifying the MSigDB category (e.g., `"H"`
#'   for Hallmark, `"C2"` for curated gene sets). Default is `"H"`.
#' @param subcategory Character string specifying the MSigDB subcategory to filter
#'   gene sets. Default is `NULL`.
#' @param palette_bar Character string or integer specifying the color palette for
#'   bar plots. Default is `"jama"`.
#' @param palette_gsea Integer specifying the color palette for GSEA plots. Default
#'   is `2`.
#' @param cols_gsea Character vector specifying custom colors for GSEA enrichment
#'   plots. If `NULL`, colors are automatically generated. Default is `NULL`.
#' @param cols_bar Character vector specifying custom colors for the enrichment
#'   bar plot. If `NULL`, colors are automatically generated. Default is `NULL`.
#' @param show_bar Integer specifying the number of top enriched gene sets to display
#'   in the bar plot. Default is `10`.
#' @param show_col Logical indicating whether to display color names in the bar plot.
#'   Default is `FALSE`.
#' @param show_plot Logical indicating whether to display GSEA enrichment plots.
#'   Default is `FALSE`.
#' @param show_gsea Integer specifying the number of top significant gene sets for
#'   which to generate GSEA plots. Default is `8`.
#' @param show_path_n Integer specifying the number of pathways to display in GSEA
#'   plots. Default is `20`.
#' @param plot_single_sig Logical indicating whether to generate separate plots for
#'   each significant gene set. Default is `TRUE`.
#' @param project Character string specifying the project name for output files.
#'   Default is `"custom_sig"`.
#' @param minGSSize Integer specifying the minimum gene set size for analysis.
#'   Default is `10`.
#' @param maxGSSize Integer specifying the maximum gene set size for analysis.
#'   Default is `500`.
#' @param verbose Logical indicating whether to display progress messages. Default
#'   is `TRUE`.
#' @param seed Logical indicating whether to set a random seed for reproducibility.
#'   Default is `FALSE`.
#' @param fig.type Character string specifying the file format for saving plots
#'   (e.g., `"pdf"`, `"png"`). Default is `"pdf"`.
#' @param print_bar Logical indicating whether to save and print the bar plot.
#'   Default is `TRUE`.
#'
#' @return List containing:
#' \describe{
#'   \item{up}{Data frame of up-regulated enriched gene sets}
#'   \item{down}{Data frame of down-regulated enriched gene sets}
#'   \item{all}{Complete GSEA results}
#'   \item{plot_top}{GSEA enrichment plot for top gene sets}
#' }
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' # Load example data
#' eset_stad <- load_data("eset_stad")
#' stad_group <- load_data("stad_group")
#' oldwd <- getwd()
#' on.exit(setwd(oldwd))
#' setwd(tempdir())
#' deg <- iobr_deg(
#'   eset = eset_stad, pdata = stad_group, group_id = "subtype",
#'   pdata_id = "ID", array = FALSE, method = "DESeq2",
#'   contrast = c("EBV", "GS"), path = "STAD"
#' )
#' # Run GSEA with custom gene sets
#' signature_tme <- load_data("signature_tme")
#' res <- sig_gsea(deg = deg, genesets = signature_tme)
#' }
sig_gsea <- function(deg,
                     genesets = NULL,
                     path = NULL,
                     gene_symbol = "symbol",
                     logfc = "log2FoldChange",
                     org = c("hsa", "mus"),
                     msigdb = TRUE,
                     category = "H",
                     subcategory = NULL,
                     palette_bar = "jama",
                     palette_gsea = 2,
                     cols_gsea = NULL,
                     cols_bar = NULL,
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
  org <- rlang::arg_match(org)

  # Input validation
  if (is.null(deg) || !is.data.frame(deg)) {
    cli::cli_abort("{.arg deg} must be a data frame with differential expression results")
  }
  if (!gene_symbol %in% colnames(deg)) {
    cli::cli_abort("Gene symbol column {.val {gene_symbol}} not found in deg")
  }
  if (!logfc %in% colnames(deg)) {
    cli::cli_abort("Log fold change column {.val {logfc}} not found in deg")
  }
  if (minGSSize < 1 || maxGSSize < minGSSize) {
    cli::cli_abort("Invalid gene set size parameters")
  }

  rlang::check_installed("clusterProfiler")

  # Setup output path
  save_results <- !is.null(path)
  abspath <- NULL
  file_store <- NULL
  if (save_results) {
    file_store <- path
    if (!dir.exists(file_store)) dir.create(file_store, recursive = TRUE)
    abspath <- file.path(normalizePath(file_store, winslash = "/", mustWork = FALSE), "")
  }

  # Safely rename columns
  if (gene_symbol != "symbol" && "symbol" %in% colnames(deg)) {
    colnames(deg)[colnames(deg) == "symbol"] <- "symbol_subs"
  }
  colnames(deg)[colnames(deg) == gene_symbol] <- "symbol"
  colnames(deg)[colnames(deg) == logfc] <- "logfc"

  # Map gene symbols to Entrez IDs
  cli::cli_alert_info("Mapping gene symbols to Entrez IDs for {.val {org}}")
  database <- switch(org,
    hsa = "org.Hs.eg.db",
    mus = "org.Mm.eg.db"
  )
  rlang::check_installed(database)

  entrizid <- clusterProfiler::bitr(deg$symbol,
    fromType = "SYMBOL",
    toType = c("GENENAME", "ENTREZID"),
    OrgDb = database
  )

  # Merge and prepare genelist
  deg <- merge(deg, entrizid, by.x = "symbol", by.y = "SYMBOL", all.x = TRUE, all.y = FALSE)
  gene_id_logfc <- deg %>%
    dplyr::select("ENTREZID", "logfc") %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(ENTREZID), !is.na(logfc)) %>%
    dplyr::arrange(dplyr::desc(logfc))

  genelist <- gene_id_logfc$logfc
  names(genelist) <- gene_id_logfc$ENTREZID
  genelist <- genelist[order(genelist, decreasing = TRUE)]

  # Get gene sets
  term2genes <- .get_gene_sets(genesets, org, category, subcategory, project)

  # Run GSEA
  cli::cli_alert_info("Running GSEA analysis...")
  hall_gsea <- clusterProfiler::GSEA(
    genelist,
    exponent = 1,
    eps = 0,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    TERM2GENE = term2genes,
    TERM2NAME = NA,
    verbose = verbose,
    seed = seed,
    by = "fgsea"
  )

  # Convert to readable format
  rlang::check_installed("DOSE")
  hall_gsea <- DOSE::setReadable(hall_gsea, OrgDb = database, keyType = "ENTREZID")

  # Save results
  if (save_results) {
    csv_file <- paste0(abspath, "1-", category, "_GSEA_significant_results.csv")
    utils::write.csv(as.data.frame(hall_gsea), file = csv_file, row.names = FALSE)
    cli::cli_alert_success("GSEA results written to: {.path {csv_file}}")
  }

  # Generate plots if results exist
  hall_gsea_result <- .process_gsea_results(
    hall_gsea, category, show_gsea, show_path_n, plot_single_sig,
    palette_gsea, cols_gsea, palette_bar, cols_bar, show_plot,
    print_bar, save_results, file_store, fig.type
  )

  if (save_results) {
    save(hall_gsea_result, file = paste0(abspath, "0-", category, "_combined_GSEA_result.RData"))
  }

  hall_gsea_result
}

#' Get gene sets for GSEA
#' @keywords internal
#' @noRd
.get_gene_sets <- function(genesets, org, category, subcategory, project) {
  if (is.null(genesets)) {
    rlang::check_installed("msigdbr")
    species <- switch(org,
      hsa = "Homo sapiens",
      mus = "Mus musculus"
    )

    m_df <- msigdbr::msigdbr(species = species)
    a <- m_df %>%
      dplyr::distinct(.data$gs_cat, .data$gs_subcat) %>%
      dplyr::arrange(.data$gs_cat, .data$gs_subcat)
    cli::cli_alert_info("Available MSigDB categories:")
    print(as.data.frame(a))

    category <- category %||% "H"
    term2genes <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
    term2genes %>%
      dplyr::select("gs_name", "entrez_gene") %>%
      as.data.frame()
  } else {
    # Process custom gene sets
    term2genes <- output_sig(genesets, file.name = "sig")
    file.remove("sig.csv")
    term2genes <- reshape2::melt(as.matrix(term2genes))
    term2genes <- term2genes[, -1]
    colnames(term2genes) <- c("gs_name", "symbol")

    database <- switch(org,
      hsa = "org.Hs.eg.db",
      mus = "org.Mm.eg.db"
    )
    rlang::check_installed(database)

    entrizid <- clusterProfiler::bitr(term2genes$symbol,
      fromType = "SYMBOL",
      toType = c("GENENAME", "ENTREZID"),
      OrgDb = database
    )

    term2genes <- merge(term2genes, entrizid, by.x = "symbol", by.y = "SYMBOL", all.x = TRUE, all.y = FALSE)
    term2genes <- term2genes[, c("gs_name", "ENTREZID")]
    colnames(term2genes) <- c("gs_name", "entrez_gene")
    term2genes
  }
}

#' Process GSEA results and generate plots
#' @keywords internal
#' @noRd
.process_gsea_results <- function(hall_gsea, category, show_gsea, show_path_n,
                                  plot_single_sig, palette_gsea, cols_gsea,
                                  palette_bar, cols_bar, show_plot, print_bar,
                                  save_results, file_store, fig.type) {
  if (nrow(hall_gsea) == 0) {
    cli::cli_alert_warning("No terms enriched under pvalue cutoff = 0.05")
    return(hall_gsea)
  }

  # Select top pathways
  n_paths <- min(show_gsea, nrow(hall_gsea))
  paths <- rownames(hall_gsea@result)[1:n_paths]
  cli::cli_alert_info("Most significant gene sets: {paste(paths, collapse = ', ')}")

  # Get colors
  gseacol <- cols_gsea %||% get_cols(palette = palette_gsea, show_col = FALSE)

  # Generate top GSEA plot
  rlang::check_installed("enrichplot")
  GSEAPLOT <- enrichplot::gseaplot2(
    x = hall_gsea,
    geneSetID = paths,
    subplots = c(1, 2),
    color = gseacol[seq_along(paths)],
    pvalue_table = TRUE
  )

  if (save_results) {
    ggplot2::ggsave(
      filename = paste0("2-", category, "_Top_", show_gsea, "_GSEA_plot.", fig.type),
      plot = GSEAPLOT, path = file_store,
      width = 11, height = 7, dpi = 300
    )
  }
  if (show_plot) print(GSEAPLOT)

  # Generate individual plots
  if (plot_single_sig) {
    paths_all <- rownames(hall_gsea@result)
    paths2 <- paths_all[1:min(show_path_n, length(paths_all))]
    paths2 <- paths2[!is.na(paths2)]

    for (i in seq_along(paths2)) {
      single_path <- paths2[i]
      cli::cli_alert_info("Processing: {.val {single_path}}")

      gseaplot <- enrichplot::gseaplot2(
        x = hall_gsea,
        geneSetID = single_path,
        subplots = c(1, 2),
        color = gseacol[i],
        pvalue_table = TRUE
      )

      if (show_plot) print(gseaplot)
      if (save_results) {
        ggplot2::ggsave(
          filename = paste0(i + 4, "-GSEA_plot-", single_path, ".", fig.type),
          plot = gseaplot, path = file_store,
          width = 11, height = 7.5, dpi = 300
        )
      }
    }
  }

  # Split into up and down
  down_gogo <- hall_gsea[hall_gsea$pvalue < 0.05 & hall_gsea$enrichmentScore < 0, ]
  if (nrow(down_gogo) > 0) down_gogo$group <- -1
  down_gogo <- down_gogo[1:min(show_path_n, nrow(down_gogo)), ]
  down_gogo <- down_gogo[!is.na(down_gogo$p.adjust), ]

  up_gogo <- hall_gsea[hall_gsea$pvalue < 0.05 & hall_gsea$enrichmentScore > 0, ]
  if (nrow(up_gogo) > 0) up_gogo$group <- 1
  up_gogo <- up_gogo[1:min(show_path_n, nrow(up_gogo)), ]
  up_gogo <- up_gogo[!is.na(up_gogo$p.adjust), ]

  # Generate bar plot
  if (print_bar) {
    gsea_bar <- enrichment_barplot(
      up_terms = up_gogo,
      down_terms = down_gogo,
      palette = palette_bar,
      cols = cols_bar,
      title = "GSEA-Enrichment",
      width_wrap = 30
    )

    n_bar <- nrow(up_gogo) + nrow(down_gogo)
    height_bar <- 0.5 * n_bar + 3

    if (save_results) {
      ggplot2::ggsave(
        filename = paste0("3-", category, "_GSEA_barplot.", fig.type),
        plot = gsea_bar, path = file_store,
        width = 6, height = height_bar
      )
    }
    if (show_plot) print(gsea_bar)
  }

  list(
    up = hall_gsea[hall_gsea$pvalue < 0.05 & hall_gsea$enrichmentScore > 0, ],
    down = hall_gsea[hall_gsea$pvalue < 0.05 & hall_gsea$enrichmentScore < 0, ],
    all = tibble::as_tibble(hall_gsea),
    plot_top = GSEAPLOT
  )
}
