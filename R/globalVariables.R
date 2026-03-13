# 文件：R/zzz-globalVariables.R
# 作用：静默 R CMD check 的 “no visible binding for global variable” NOTE
# 用法：直接放在 R/ 目录下，无需额外 source，也无需 @export
utils::globalVariables(c(
  # ---------- 管道 / data.table / tidyverse 内部变量 ----------
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
  "FPR", "TPR", "variable","itor",

  # ---------- 跨包函数 / S3 方法 ----------
  # 以下名字在代码里通过 :: 或 get() 调用，但仍被 codetools 报 undefined
  "DefaultAssay<-", "check_installed", "geom_hdr_lines",
  "ggboxplot", "grid.arrange", "rlm", "lsei", "DESeq", "results",
  "counts", "makeContrasts", "eBayes", "lmFit", "contrasts.fit",
  "topTable", "stat_compare_means", "compare_means",
  "listDatasets", "useMart", "getLDS", "random_strata_cells",
  "bicor", "fundamentalNetworkConcepts", "anno_barplot", "gene_count"
))
