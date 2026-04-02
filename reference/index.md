# Package index

## Data Input & Validation

Read, validate, and format expression data.

- [`check_eset()`](https://iobr.github.io/IOBR/reference/check_eset.md)
  : Check Integrity and Outliers of Expression Set
- [`anno_eset()`](https://iobr.github.io/IOBR/reference/anno_eset.md) :
  Annotate Gene Expression Matrix and Remove Duplicated Genes
- [`combine_pd_eset()`](https://iobr.github.io/IOBR/reference/combine_pd_eset.md)
  : Combine Phenotype Data and Expression Set
- [`merge_eset()`](https://iobr.github.io/IOBR/reference/merge_eset.md)
  : Merge Expression Sets by Row Names
- [`rbind_iobr()`](https://iobr.github.io/IOBR/reference/rbind_iobr.md)
  : Row Bind Multiple Data Sets
- [`filterCommonGenes()`](https://iobr.github.io/IOBR/reference/filterCommonGenes.md)
  : filterCommonGenes
- [`remove_duplicate_genes()`](https://iobr.github.io/IOBR/reference/remove_duplicate_genes.md)
  : Remove Duplicate Gene Symbols in Gene Expression Data
- [`merge_duplicate()`](https://iobr.github.io/IOBR/reference/merge_duplicate.md)
  : Merge Data Frames with Duplicated Column Names
- [`assimilate_data()`](https://iobr.github.io/IOBR/reference/assimilate_data.md)
  : Harmonize Two Data Frames by Column Structure

## Data Preprocessing & Transformation

Normalize, transform, and batch-correct expression data.

- [`transform_data()`](https://iobr.github.io/IOBR/reference/transform_data.md)
  : Transform NA, Inf, or Zero Values in Data
- [`count2tpm()`](https://iobr.github.io/IOBR/reference/count2tpm.md) :
  Convert Read Counts to Transcripts Per Million (TPM)
- [`log2eset()`](https://iobr.github.io/IOBR/reference/log2eset.md) :
  Log2 Transformation of Gene Expression Matrix
- [`remove_batcheffect()`](https://iobr.github.io/IOBR/reference/remove_batcheffect.md)
  : Removing Batch Effect from Expression Sets
- [`RemoveBatchEffect()`](https://iobr.github.io/IOBR/reference/RemoveBatchEffect.md)
  : Remove Batch Effect of Expression Set
- [`scale_matrix()`](https://iobr.github.io/IOBR/reference/scale_matrix.md)
  : Scale and Manipulate a Matrix
- [`eset_distribution()`](https://iobr.github.io/IOBR/reference/eset_distribution.md)
  : Visualize Expression Set Distribution
- [`mouse2human_eset()`](https://iobr.github.io/IOBR/reference/mouse2human_eset.md)
  : Convert Mouse Gene Symbols to Human Gene Symbols
- [`patterns_to_na`](https://iobr.github.io/IOBR/reference/patterns_to_na.md)
  : Default Pattern List for Name Cleaning
- [`tcga_rna_preps()`](https://iobr.github.io/IOBR/reference/tcga_rna_preps.md)
  : Preprocess TCGA RNA-seq Data

## TME Deconvolution — Main Interface

Estimate immune and stromal cell fractions from bulk expression data
using 11 integrated deconvolution algorithms.

- [`deconvo_tme()`](https://iobr.github.io/IOBR/reference/deconvo_tme.md)
  : Main TME Deconvolution Function
- [`tme_deconvolution_methods`](https://iobr.github.io/IOBR/reference/tme_deconvolution_methods.md)
  : TME Deconvolution Methods
- [`select_method()`](https://iobr.github.io/IOBR/reference/select_method.md)
  : Select a Signature Scoring Method Subset
- [`iobr_deconvo_pipeline()`](https://iobr.github.io/IOBR/reference/iobr_deconvo_pipeline.md)
  : Tumor Microenvironment (TME) Deconvolution Pipeline

## TME Deconvolution — Individual Methods

Standalone wrappers for each supported deconvolution method.

- [`deconvo_cibersort()`](https://iobr.github.io/IOBR/reference/deconvo_cibersort.md)
  : Deconvolve Using CIBERSORT
- [`deconvo_timer()`](https://iobr.github.io/IOBR/reference/deconvo_timer.md)
  : Deconvolve Using TIMER
- [`deconvo_xcell()`](https://iobr.github.io/IOBR/reference/deconvo_xcell.md)
  : Deconvolve Immune Microenvironment Using xCell
- [`deconvo_mcpcounter()`](https://iobr.github.io/IOBR/reference/deconvo_mcpcounter.md)
  : Deconvolve Immune Microenvironment Using MCP-Counter
- [`deconvo_estimate()`](https://iobr.github.io/IOBR/reference/deconvo_estimate.md)
  : Calculate ESTIMATE Scores
- [`deconvo_epic()`](https://iobr.github.io/IOBR/reference/deconvo_epic.md)
  : Deconvolve Immune Microenvironment Using EPIC
- [`deconvo_ips()`](https://iobr.github.io/IOBR/reference/deconvo_ips.md)
  : Calculate Immunophenoscore (IPS)
- [`deconvo_quantiseq()`](https://iobr.github.io/IOBR/reference/deconvo_quantiseq.md)
  : Deconvolve Using quanTIseq
- [`deconvolute_quantiseq.default()`](https://iobr.github.io/IOBR/reference/deconvolute_quantiseq.default.md)
  : Use quanTIseq to Deconvolute a Gene Expression Matrix
- [`deconvolute_timer.default()`](https://iobr.github.io/IOBR/reference/deconvolute_timer.default.md)
  : Deconvolute Tumor Microenvironment Using TIMER
- [`deconvo_ref()`](https://iobr.github.io/IOBR/reference/deconvo_ref.md)
  : Deconvolve Using Custom Reference

## TME Deconvolution — CIBERSORT Internals

Core CIBERSORT algorithm implementation.

- [`CIBERSORT()`](https://iobr.github.io/IOBR/reference/CIBERSORT.md) :
  CIBERSORT Deconvolution Algorithm
- [`CoreAlg()`](https://iobr.github.io/IOBR/reference/CoreAlg.md) : Core
  Algorithm for CIBERSORT Deconvolution
- [`GetFractions.Abbas()`](https://iobr.github.io/IOBR/reference/GetFractions.Abbas.md)
  : Constrained Regression Method (Abbas et al., 2009)
- [`ParseInputExpression()`](https://iobr.github.io/IOBR/reference/ParseInputExpression.md)
  : Parse Input Gene Expression Data
- [`doPerm()`](https://iobr.github.io/IOBR/reference/doPerm.md) :
  Permutation Test for CIBERSORT
- [`parallel_doperm()`](https://iobr.github.io/IOBR/reference/parallel_doperm.md)
  : Parallel Permutation Test for CIBERSORT

## TME Deconvolution — TIMER Internals

TIMER algorithm utilities and cancer-type support.

- [`timer_available_cancers`](https://iobr.github.io/IOBR/reference/timer_available_cancers.md)
  : TIMER Available Cancer Types
- [`timer_info()`](https://iobr.github.io/IOBR/reference/timer_info.md)
  : Source code for the TIMER deconvolution method.
- [`check_cancer_types()`](https://iobr.github.io/IOBR/reference/check_cancer_types.md)
  : Process Batch Table and Validate Cancer Types
- [`GetOutlierGenes()`](https://iobr.github.io/IOBR/reference/GetOutlierGenes.md)
  : Get Outlier Genes
- [`Top_probe()`](https://iobr.github.io/IOBR/reference/Top_probe.md) :
  Top Probe Selector

## TME Deconvolution — Other Algorithm Internals

Internal helpers for EPIC, IPS, MCPcounter, and ESTIMATE.

- [`MCPcounter.estimate()`](https://iobr.github.io/IOBR/reference/MCPcounter.estimate.md)
  : MCP-counter Cell Population Abundance Estimation
- [`EPIC()`](https://iobr.github.io/IOBR/reference/EPIC.md) : Estimate
  the proportion of immune and cancer cells.
- [`IPS_calculation()`](https://iobr.github.io/IOBR/reference/IPS_calculation.md)
  : Calculate Immunophenoscore (IPS)
- [`ipsmap()`](https://iobr.github.io/IOBR/reference/ipsmap.md) : Map
  Score to Immunophenoscore
- [`estimateScore()`](https://iobr.github.io/IOBR/reference/estimateScore.md)
  : estimateScore
- [`plotPurity()`](https://iobr.github.io/IOBR/reference/plotPurity.md)
  : plotPurity

## Signature Scoring

Calculate gene signature scores using PCA, z-score, and ssGSEA methods
for 300+ curated TME signatures.

- [`calculate_sig_score()`](https://iobr.github.io/IOBR/reference/calculate_sig_score.md)
  : Calculate Signature Score
- [`sigScore()`](https://iobr.github.io/IOBR/reference/sigScore.md) :
  Calculate Signature Score Using PCA, Mean, or Z-score Methods
- [`signature_score_calculation_methods`](https://iobr.github.io/IOBR/reference/signature_score_calculation_methods.md)
  : Signature Score Calculation Methods
- [`calculate_sig_score_pca()`](https://iobr.github.io/IOBR/reference/calculate_sig_score_pca.md)
  : Calculate Signature Score Using PCA Method
- [`calculate_sig_score_zscore()`](https://iobr.github.io/IOBR/reference/calculate_sig_score_zscore.md)
  : Calculate Signature Score Using Z-Score Method
- [`calculate_sig_score_ssgsea()`](https://iobr.github.io/IOBR/reference/calculate_sig_score_ssgsea.md)
  : Calculate Signature Score Using ssGSEA Method
- [`calculate_sig_score_integration()`](https://iobr.github.io/IOBR/reference/calculate_sig_score_integration.md)
  : Calculate Signature Score Using Integration Method
- [`test_for_infiltration()`](https://iobr.github.io/IOBR/reference/test_for_infiltration.md)
  : Test for Cell Population Infiltration

## Statistical Analysis — Correlation

Batch correlation tests between features or signature scores.

- [`batch_cor()`](https://iobr.github.io/IOBR/reference/batch_cor.md) :
  Batch Correlation Analysis
- [`batch_pcc()`](https://iobr.github.io/IOBR/reference/batch_pcc.md) :
  Batch Calculation of Partial Correlation Coefficients
- [`get_cor()`](https://iobr.github.io/IOBR/reference/get_cor.md) :
  Calculate and Visualize Correlation Between Two Variables
- [`get_cor_matrix()`](https://iobr.github.io/IOBR/reference/get_cor_matrix.md)
  : Calculate and Visualize Correlation Matrix Between Two Variable Sets

## Statistical Analysis — Group Comparison

Statistical tests comparing groups (Wilcoxon, Kruskal-Wallis).

- [`batch_wilcoxon()`](https://iobr.github.io/IOBR/reference/batch_wilcoxon.md)
  : Batch Wilcoxon Rank-Sum Test Between Two Groups
- [`batch_kruskal()`](https://iobr.github.io/IOBR/reference/batch_kruskal.md)
  : Batch Kruskal-Wallis Test
- [`exact_pvalue()`](https://iobr.github.io/IOBR/reference/exact_pvalue.md)
  : Calculate Exact P-Value for Correlation

## Statistical Analysis — Survival & Time-to-Event

Optimal cutoff identification and time-dependent ROC analysis.

- [`best_cutoff()`](https://iobr.github.io/IOBR/reference/best_cutoff.md)
  : Extract Best Cutoff and Add Binary Variable to Data Frame
- [`best_cutoff2()`](https://iobr.github.io/IOBR/reference/best_cutoff2.md)
  : Extract Best Cutoff and Add Binary Variable to Data Frame
- [`calculate_break_month()`](https://iobr.github.io/IOBR/reference/calculate_break_month.md)
  : Break Time Into Blocks
- [`CalculateTimeROC()`](https://iobr.github.io/IOBR/reference/CalculateTimeROC.md)
  : Calculate Time-Dependent ROC Curve
- [`PlotTimeROC()`](https://iobr.github.io/IOBR/reference/PlotTimeROC.md)
  : Plot Time-Dependent ROC Curves
- [`roc_time()`](https://iobr.github.io/IOBR/reference/roc_time.md) :
  Time-dependent ROC Curve for Survival Analysis

## Visualization — Heatmaps & Correlation

Heatmaps and correlation plots for TME features.

- [`sig_heatmap()`](https://iobr.github.io/IOBR/reference/sig_heatmap.md)
  : Signature Heatmap with Optional Annotations
- [`sig_pheatmap()`](https://iobr.github.io/IOBR/reference/sig_pheatmap.md)
  : Generate Heatmap for Signature Data
- [`iobr_cor_plot()`](https://iobr.github.io/IOBR/reference/iobr_cor_plot.md)
  : Integrative Correlation Analysis Between Phenotype and Features

## Visualization — Survival & ROC

Kaplan-Meier survival curves and ROC plots.

- [`sig_surv_plot()`](https://iobr.github.io/IOBR/reference/sig_surv_plot.md)
  : Generate Kaplan-Meier Survival Plot for Signature
- [`batch_surv()`](https://iobr.github.io/IOBR/reference/batch_surv.md)
  : Batch Survival Analysis
- [`batch_sig_surv_plot()`](https://iobr.github.io/IOBR/reference/batch_sig_surv_plot.md)
  : Batch Signature Survival Plot
- [`surv_group()`](https://iobr.github.io/IOBR/reference/surv_group.md)
  : Generate Kaplan-Meier Survival Plots for Categorical Groups
- [`subgroup_survival()`](https://iobr.github.io/IOBR/reference/subgroup_survival.md)
  : Subgroup Survival Analysis Using Cox Proportional Hazards Models
- [`sig_roc()`](https://iobr.github.io/IOBR/reference/sig_roc.md) : Plot
  ROC Curves and Compare Them

## Visualization — Comparison Plots

Box plots, violin plots, and forest plots.

- [`sig_box()`](https://iobr.github.io/IOBR/reference/sig_box.md) :
  Signature Box Plot with Statistical Comparisons
- [`sig_box_batch()`](https://iobr.github.io/IOBR/reference/sig_box_batch.md)
  : Batch Signature Box Plots for Group Comparisons
- [`sig_forest()`](https://iobr.github.io/IOBR/reference/sig_forest.md)
  : Forest Plot for Survival Analysis Results

## Visualization — Distribution & Enrichment

Bar charts, pie charts, and enrichment visualizations.

- [`cell_bar_plot()`](https://iobr.github.io/IOBR/reference/cell_bar_plot.md)
  : Visualize Cell Fractions as Stacked Bar Chart
- [`percent_bar_plot()`](https://iobr.github.io/IOBR/reference/percent_bar_plot.md)
  : Create a Percent Bar Plot
- [`pie_chart()`](https://iobr.github.io/IOBR/reference/pie_chart.md) :
  Create Pie or Donut Charts
- [`enrichment_barplot()`](https://iobr.github.io/IOBR/reference/enrichment_barplot.md)
  : Enrichment Bar Plot with Two Directions

## Visualization — Themes & Color Utilities

Color palettes, themes, and plot formatting helpers.

- [`get_cols()`](https://iobr.github.io/IOBR/reference/get_cols.md) :
  Set and View Color Palettes
- [`mapcolors()`](https://iobr.github.io/IOBR/reference/mapcolors.md) :
  Map Score to Color
- [`mapbw()`](https://iobr.github.io/IOBR/reference/mapbw.md) : Map
  Score to Black and White Color
- [`design_mytheme()`](https://iobr.github.io/IOBR/reference/design_mytheme.md)
  : Design Custom Theme for ggplot2 Plots
- [`palettes()`](https://iobr.github.io/IOBR/reference/palettes.md) :
  Select Color Palettes for Visualization
- [`DrawQQPlot()`](https://iobr.github.io/IOBR/reference/DrawQQPlot.md)
  : Draw QQ Plot Comparing Cancer and Immune Expression

## Prognostic Modeling

Build LASSO/Elastic-Net prognostic models and compute risk scores.

- [`PrognosticModel()`](https://iobr.github.io/IOBR/reference/PrognosticModel.md)
  : Build Prognostic Models Using LASSO and Ridge Regression
- [`BinomialModel()`](https://iobr.github.io/IOBR/reference/BinomialModel.md)
  : Binomial Model Construction
- [`add_riskscore()`](https://iobr.github.io/IOBR/reference/add_riskscore.md)
  : Add Risk Score to Dataset
- [`PrognosticAUC()`](https://iobr.github.io/IOBR/reference/PrognosticAUC.md)
  : Calculate Time-Dependent AUC for Survival Models
- [`BinomialAUC()`](https://iobr.github.io/IOBR/reference/BinomialAUC.md)
  : Calculate AUC for Binomial Model
- [`SplitTrainTest()`](https://iobr.github.io/IOBR/reference/SplitTrainTest.md)
  : Split Data into Training and Testing Sets
- [`PrognosticResult()`](https://iobr.github.io/IOBR/reference/PrognosticResult.md)
  : Compute Prognostic Results for Survival Models
- [`RegressionResult()`](https://iobr.github.io/IOBR/reference/RegressionResult.md)
  : Regression Result Computation
- [`PlotAUC()`](https://iobr.github.io/IOBR/reference/PlotAUC.md) : Plot
  AUC ROC Curves
- [`getHRandCIfromCoxph()`](https://iobr.github.io/IOBR/reference/getHRandCIfromCoxph.md)
  : Extract Hazard Ratio and Confidence Intervals from Cox Model
- [`ProcessingData()`](https://iobr.github.io/IOBR/reference/ProcessingData.md)
  : Process Data for Model Construction
- [`Enet()`](https://iobr.github.io/IOBR/reference/Enet.md) : Elastic
  Net Model Fitting
- [`CalculatePref()`](https://iobr.github.io/IOBR/reference/CalculatePref.md)
  : Calculate Performance Metrics

## Clustering & Subgroup Analysis

Cluster TME profiles and identify patient subgroups.

- [`tme_cluster()`](https://iobr.github.io/IOBR/reference/tme_cluster.md)
  : Identification of TME Cluster

## Feature Selection

Identify informative features and variable genes.

- [`feature_select()`](https://iobr.github.io/IOBR/reference/feature_select.md)
  : Feature Selection via Correlation or Differential Expression
- [`feature_manipulation()`](https://iobr.github.io/IOBR/reference/feature_manipulation.md)
  : Feature Quality Control and Filtering
- [`lasso_select()`](https://iobr.github.io/IOBR/reference/lasso_select.md)
  : Feature Selection for Predictive or Prognostic Models Using LASSO
  Regression
- [`high_var_fea()`](https://iobr.github.io/IOBR/reference/high_var_fea.md)
  : Identify High-Variance Features from Statistical Results
- [`find_variable_genes()`](https://iobr.github.io/IOBR/reference/find_variable_genes.md)
  : Identify Variable Genes in Expression Data
- [`find_outlier_samples()`](https://iobr.github.io/IOBR/reference/find_outlier_samples.md)
  : Identify Outlier Samples in Gene Expression Data

## Genomics & Differential Expression

Differential expression, mutation analysis, and PCA.

- [`limma.dif()`](https://iobr.github.io/IOBR/reference/limma.dif.md) :
  Differential Expression Analysis Using Limma
- [`iobr_deg()`](https://iobr.github.io/IOBR/reference/iobr_deg.md) :
  Differential Expression Analysis
- [`iobr_pca()`](https://iobr.github.io/IOBR/reference/iobr_pca.md) :
  Principal Component Analysis (PCA) Visualization
- [`find_mutations()`](https://iobr.github.io/IOBR/reference/find_mutations.md)
  : Analyze Mutations Related to Signature Scores
- [`find_markers_in_bulk()`](https://iobr.github.io/IOBR/reference/find_markers_in_bulk.md)
  : Identify Marker Features in Bulk Expression Data
- [`make_mut_matrix()`](https://iobr.github.io/IOBR/reference/make_mut_matrix.md)
  : Construct Mutation Matrices from MAF Data
- [`ConvertRownameToLoci()`](https://iobr.github.io/IOBR/reference/ConvertRownameToLoci.md)
  : Convert Rowname To Loci

## scRNA-seq Reference Construction

Build custom deconvolution reference matrices from scRNA-seq data.

- [`generateRef()`](https://iobr.github.io/IOBR/reference/generateRef.md)
  : Generate Reference Signature Matrix
- [`generateRef_rnaseq()`](https://iobr.github.io/IOBR/reference/generateRef_rnaseq.md)
  : Generate Reference Gene Matrix from RNA-seq DEGs
- [`generateRef_seurat()`](https://iobr.github.io/IOBR/reference/generateRef_seurat.md)
  : Generate Reference Matrix from Seurat Object
- [`generateRef_limma()`](https://iobr.github.io/IOBR/reference/generateRef_limma.md)
  : Generate Reference Signature Matrix Using Limma
- [`generateRef_DEseq2()`](https://iobr.github.io/IOBR/reference/generateRef_DEseq2.md)
  : Generate Reference Signature Matrix Using DESeq2
- [`extract_sc_data()`](https://iobr.github.io/IOBR/reference/extract_sc_data.md)
  : Extract Data Frame from Seurat Object
- [`get_sig_sc()`](https://iobr.github.io/IOBR/reference/get_sig_sc.md)
  : Extract Top Marker Genes from Single-Cell Differential Results
- [`Construct_con()`](https://iobr.github.io/IOBR/reference/Construct_con.md)
  : Construct Contrast Matrix

## Signature Format & Export

Format, annotate, and export gene signatures.

- [`format_signatures()`](https://iobr.github.io/IOBR/reference/format_signatures.md)
  : Transform Signature Data into List Format
- [`format_msigdb()`](https://iobr.github.io/IOBR/reference/format_msigdb.md)
  : Format Input Signatures from MSigDB
- [`output_sig()`](https://iobr.github.io/IOBR/reference/output_sig.md)
  : Save Signature Data to File
- [`outputGCT()`](https://iobr.github.io/IOBR/reference/outputGCT.md) :
  outputGCT
- [`sig_gsea()`](https://iobr.github.io/IOBR/reference/sig_gsea.md) :
  Perform Gene Set Enrichment Analysis (GSEA)

## Ligand-Receptor Interactions

Calculate ligand-receptor interaction scores.

- [`LR_cal()`](https://iobr.github.io/IOBR/reference/LR_cal.md) :
  Calculate Ligand-Receptor Interaction Scores

## General Utilities

Miscellaneous helper functions for data and file management.

- [`load_data()`](https://iobr.github.io/IOBR/reference/load_data.md) :
  Load IOBR Datasets
- [`creat_folder()`](https://iobr.github.io/IOBR/reference/creat_folder.md)
  : Create Nested Output Folders
- [`remove_names()`](https://iobr.github.io/IOBR/reference/remove_names.md)
  : Remove Patterns from Column Names or Variables

## Datasets — Expression Sets

Example bulk RNA-seq expression matrices.

- [`eset_stad`](https://iobr.github.io/IOBR/reference/eset_stad.md) :
  Toy STAD expression matrix
- [`eset_blca`](https://iobr.github.io/IOBR/reference/eset_blca.md) :
  TCGA-BLCA Bladder Cancer Expression Data
- [`eset_gse62254`](https://iobr.github.io/IOBR/reference/eset_gse62254.md)
  : GSE62254 Gastric Cancer Expression Data
- [`eset_tme_stad`](https://iobr.github.io/IOBR/reference/eset_tme_stad.md)
  : TCGA-STAD Tumor Microenvironment Signature Scores

## Datasets — Deconvolution Reference Matrices

Reference signature matrices for CIBERSORT and EPIC.

- [`lm22`](https://iobr.github.io/IOBR/reference/lm22.md) : Reference
  profiles for immune cell types using lm22 (EPIC/IOBR)
- [`TRef`](https://iobr.github.io/IOBR/reference/TRef.md) : Reference
  profiles from tumor-infiltrating non-malignant cells
- [`BRef`](https://iobr.github.io/IOBR/reference/BRef.md) : Reference
  profiles for B cell–related deconvolution (EPIC/IOBR)

## Datasets — Gene Signatures

Curated gene signature collections for TME scoring.

- [`signature_collection`](https://iobr.github.io/IOBR/reference/signature_collection.md)
  : Gene signature collection for pathway and immune analysis
- [`imvigor210_sig`](https://iobr.github.io/IOBR/reference/imvigor210_sig.md)
  : IMvigor210 Bladder Cancer Cohort Multi-omics Signatures
- [`sig_stad`](https://iobr.github.io/IOBR/reference/sig_stad.md) :
  TCGA-STAD Gastric Cancer Cohort with Molecular and Clinical Data
- [`tcga_stad_sig`](https://iobr.github.io/IOBR/reference/tcga_stad_sig.md)
  : TCGA-STAD Gastric Cancer Immune Infiltration Signatures

## Datasets — Sample Metadata

Clinical and phenotypic data for example datasets.

- [`pdata_stad`](https://iobr.github.io/IOBR/reference/pdata_stad.md) :
  Toy STAD Phenotype Data
- [`tcga_stad_pdata`](https://iobr.github.io/IOBR/reference/tcga_stad_pdata.md)
  : TCGA-STAD Clinical and Molecular Annotation Data
- [`imvigor210_pdata`](https://iobr.github.io/IOBR/reference/imvigor210_pdata.md)
  : IMvigor210 Bladder Cancer Immunotherapy Cohort Data
- [`sig_group`](https://iobr.github.io/IOBR/reference/sig_group.md) :
  Grouped gene signatures for IOBR analysis
- [`stad_group`](https://iobr.github.io/IOBR/reference/stad_group.md) :
  Example Clinical Data for TCGA-STAD Gastric Cancer Analysis
- [`subgroup_data`](https://iobr.github.io/IOBR/reference/subgroup_data.md)
  : Example Dataset for Subgroup Survival Analysis

## Datasets — Genomic Annotation

Gene and probe annotation tables for multiple platforms.

- [`anno_rnaseq`](https://iobr.github.io/IOBR/reference/anno_rnaseq.md)
  : General RNA-seq Annotation
- [`anno_grch38`](https://iobr.github.io/IOBR/reference/anno_grch38.md)
  : GRCh38 Human Genome Annotation
- [`anno_illumina`](https://iobr.github.io/IOBR/reference/anno_illumina.md)
  : Illumina Microarray Annotation
- [`anno_hug133plus2`](https://iobr.github.io/IOBR/reference/anno_hug133plus2.md)
  : Affymetrix Human Genome U133 Plus 2.0 Array Annotation
- [`anno_gc_vm32`](https://iobr.github.io/IOBR/reference/anno_gc_vm32.md)
  : Mouse Genome Annotation (GC/VM32)

## Datasets — Miscellaneous

Additional reference data used by package functions.

- [`deg`](https://iobr.github.io/IOBR/reference/deg.md) : Single-cell
  RNA-seq Differential Expression Analysis Results
- [`null_models`](https://iobr.github.io/IOBR/reference/null_models.md)
  : NULL Model Coefficients for MCPcounter
