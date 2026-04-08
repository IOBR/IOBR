# R/globalVariables.R
# Purpose: Suppress "no visible binding for global variable" NOTES from R CMD check
# Usage: Place in R/ directory; automatically sourced, no @export needed

#' Global Variables Declaration
#'
#' @description
#' This file declares global variables used in the IOBR package to suppress
#' R CMD check NOTES about "no visible binding for global variable". These
#' variables are primarily used in tidyverse/dplyr pipelines and data.table
#' operations where non-standard evaluation (NSE) makes static analysis difficult.
#'
#' @details
#' The variables listed here fall into several categories:
#' - tidyverse NSE variables: column names used in dplyr/data.table pipelines
#' - Internal data: package datasets referenced in functions
#' - Cross-package S3 methods: generic functions from suggested packages
#'
#' @keywords internal
#' @noRd
utils::globalVariables(c(
  # ---------- tidyverse / data.table internal variables ----------
  ".", "..p.format..", "value", "x", "y", "p.value", "p.adj", "padj",
  "label", "group1", "group2", "mean_group1", "mean_group2",
  "Estimate", "logFC", "adj.P.Val", "auc", "freq", "Freq",
  "percent_weight", "lab.ypos", "Prop", "P", "cell_type", "fraction",
  "ID", "ID1", "ID2", "Intensity", "Sample.Name", "Sample.Num",
  "Z.score", "title", "target_group", "sig_group", "variables",
  "idd", "goup3", "categorys", "point.alpha", "mutation",
  "entrez_gene", "gs_cat", "gs_name", "gs_subcat", "hallmark",
  "go_bp", "go_cc", "go_mf", "kegg", "reactome",
  "signature_collection", "immuneCuratedData", "cancer_type_genes",
  "SI_geneset", "TRef", "mRNA_cell_default", "ips_gene_set", "y1", "y2",
  "lm22", "anno_grch38", "anno_gc_vm32", "mcp_genes", "mcp_probesets",
  "quantiseq_data", "xCell.data", "null_models", "common_genes",
  "avg_log2FC", "mus_human_gene_symbol", "PurityDataAffy",
  "my_palette", "my_palette2", "Group", "ENTREZID", "logfc", "p",
  "FPR", "TPR", "variable", "itor", "s",

  # ---------- Cross-package functions / S3 methods ----------
  # Names called via :: or get() but still flagged by codetools

  # ---------- Variables added during code quality improvements ----------
  "main", "p.format", ".break_month", "draw", "coef", "Mode", "tb"
))
