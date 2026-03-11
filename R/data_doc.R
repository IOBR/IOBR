#' Toy STAD expression matrix
#'
#' A subset of TCGA-STAD RNA-seq count data used in the IOBR vignette and unit
#' tests. Rows correspond to Ensembl gene IDs and columns correspond to TCGA
#' sample barcodes. Values are raw (unnormalized) read counts.
#'
#' @format A numeric matrix with genes in rows and samples in columns.
#' @usage data(eset_stad)
#' @keywords datasets
"eset_stad"


#' Toy STAD Phenotype Data
#'
#' A data frame containing clinical and pathological annotations for the
#' TCGA stomach adenocarcinoma (STAD) cohort. Each row corresponds to one
#' tumour sample and can be matched to the columns of \code{eset_stad} via
#' the \code{ID} column. This dataset is typically used together with
#' \code{eset_stad} in examples of survival analysis, subgroup comparison
#' and immune deconvolution in the IOBR package.
#'
#' @format A data frame with one row per TCGA-STAD sample and 8 variables:
#' \describe{
#'   \item{ID}{Character. TCGA sample barcode, matching the column names
#'   of \code{eset_stad}.}
#'   \item{stage}{Factor. Pathological stage (e.g. \code{"Stage_I"},
#'   \code{"Stage_II"}, \code{"Stage_III"}, \code{"Stage_IV"}).}
#'   \item{status}{Factor. Vital status at last follow-up
#'   (\code{"Alive"} or \code{"Dead"}).}
#'   \item{Lauren}{Factor. Lauren classification of gastric cancer
#'   (\code{"Intestinal"}, \code{"Diffuse"}, \code{"Mixed"} or \code{NA}).}
#'   \item{subtype}{Factor. Molecular subtype (e.g. \code{"CIN"}, \code{"EBV"},
#'   \code{"GS"}, \code{"MSI"}).}
#'   \item{EBV}{Factor. EBV status of the tumour (\code{"Positive"} or
#'   \code{"Negative"}).}
#'   \item{time}{Numeric. Overall survival or follow-up time, typically measured in months.}
#'   \item{OS_status}{Integer/binary. Overall survival status indicator. 
#'   (\code{1} = death, \code{0} = censored)}
#' }
#' 
#' @usage data(pdata_stad)
#' @keywords datasets
"pdata_stad"




#' Gene signature collection for pathway and immune analysis
#'
#' A named list of gene signatures used in the IOBR package for immune
#' deconvolution, pathway scoring, functional annotation, and tumour
#' microenvironment (TME) characterization. Each element corresponds to a
#' predefined biological signature and contains a character vector of HGNC
#' gene symbols.
#'
#' @format A named list of length 323. Each element is a character vector of
#'   gene symbols. Representative entries include:
#'   \describe{
#'     \item{CD_8_T_effector}{Markers of CD8\eqn{^{+}} effector T cells.}
#'     \item{DDR}{DNA damage response and repair genes.}
#'     \item{Immune_Checkpoint}{Immune checkpoint molecules.}
#'     \item{CellCycle_Reg}{Core regulators of cell-cycle progression.}
#'     \item{Mismatch_Repair}{Mismatch-repair pathway genes.}
#'     \item{TMEsocreA_CIR}{TME-related signature used in TMEscore analysis.}
#'     \item{...}{Additional signatures are included in the list but are not
#'       individually listed here; all follow the same structure.}
#'   }
#'
#' @usage data(signature_collection)
#' @keywords datasets
"signature_collection"





#' Grouped gene signatures for IOBR analysis
#'
#' A named list that organizes gene signatures into functional or biological
#' categories. Each element of the list is a character vector containing the
#' names of gene signatures defined in \code{signature_collection}. A total of
#' 43 signature groups are included, covering tumour intrinsic pathways,
#' immune-related processes, stromal activity, TME characteristics and
#' immuno-oncology biomarkers. These groups are used in IOBR to conveniently
#' select sets of signatures for scoring and visualization.
#'
#' @format A named list of length 43. Each element is a character vector of
#'   signature names. Representative groups include:
#'   \describe{
#'     \item{tumor_signature}{Signatures related to intrinsic tumour biology
#'       such as cell cycle, DNA damage repair and histone regulation.}
#'     \item{EMT}{Epithelial–mesenchymal transition (EMT)–associated signatures.}
#'     \item{io_biomarkers}{Immuno-oncology biomarker–related signatures.}
#'     \item{immu_microenvironment}{Immune microenvironment–related signatures.}
#'     \item{immu_suppression}{Immune suppression–related signatures.}
#'     \item{immu_exclusion}{Signatures associated with immune exclusion and
#'       stromal barriers.}
#'     \item{TCR_BCR}{T-cell and B-cell receptor pathway signatures.}
#'     \item{tme_signatures1}{Tumour microenvironment signature panel (set 1).}
#'     \item{tme_signatures2}{Tumour microenvironment signature panel (set 2).}
#'     \item{...}{Additional groups are included (43 total), but not listed
#'       individually here; all groups follow the same structure.}
#'   }
#'
#' @usage data(sig_group)
#' @keywords datasets
"sig_group"




#' Reference profiles from tumor-infiltrating non-malignant cells
#'
#' A dataset providing reference gene expression profiles of main non-malignant
#' cell types infiltrating solid tumours, to be used with the EPIC deconvolution
#' method (and wrapped by IOBR). Expression values are based on single-cell
#' RNA-seq data from melanoma tumour samples, and are given in TPM units.
#'
#' @format A list with three components:
#' \describe{
#'   \item{refProfiles}{Numeric matrix (\code{nGenes} × \code{nCellTypes})
#'   containing the reference gene expression (TPM counts) for each gene in
#'   each cell type. Row names are gene identifiers, column names are cell
#'   types (e.g. B cells, CAFs, CD4 T, CD8 T, endothelial cells, macrophages,
#'   NK cells).}
#'
#'   \item{refProfiles.var}{Numeric matrix with the same dimensions as
#'   \code{refProfiles}, giving the variability (e.g. variance) of gene
#'   expression for each gene–cell-type combination. This matrix is used as
#'   weights in the EPIC optimisation.}
#'
#'   \item{sigGenes}{Character vector listing the signature genes used by EPIC
#'   to deconvolve cell-type proportions. Only these genes are used in the fit,
#'   although \code{refProfiles} usually stores all available genes.}
#' }
#'
#' @source Single-cell RNA-seq data of tumour-infiltrating non-malignant cells from melanoma (GEO accession \code{GSE72056}, 9 donors).  
#' @references
#' Tirosh et al. (2016) Science 352:189-196. doi:10.1126/science.aad0501
#' Racle et al. (2017) eLife 6:e26476. doi:10.7554/eLife.26476
#'
#' @usage data(TRef)
#' @keywords datasets
"TRef"






#' Reference profiles for B cell–related deconvolution (EPIC/IOBR)
#'
#' \code{BRef} provides the B-cell–related reference gene expression profiles
#' used by the EPIC deconvolution model. Within the IOBR workflow, this dataset
#' is used to estimate immune cell fractions in bulk RNA-seq samples.
#' The dataset contains normalized gene expression profiles, expression
#' variability, and the signature genes required by EPIC.
#'
#' @details
#' This reference dataset includes the following components:
#'
#' \describe{
#'   \item{\code{refProfiles}}{
#'     A numeric matrix of size \code{nGenes × nCellTypes}.
#'     Rows correspond to genes, and columns correspond to B-cell–related
#'     subtypes. The values generally represent normalized expression
#'     (e.g., TPM). EPIC uses this matrix as the baseline expression profile
#'     for the reference cell types during deconvolution.
#'   }
#'
#'   \item{\code{refProfiles.var}}{
#'     A numeric matrix with identical dimensions as \code{refProfiles}.
#'     It represents gene-level expression variability for each cell type.
#'     These variability estimates are used as weights in EPIC’s
#'     weighted least-squares optimization. If variability is not used in a
#'     given workflow, EPIC assumes the same variability for all genes.
#'   }
#'
#'   \item{\code{sigGenes}}{
#'     A character vector listing the signature genes used by EPIC for
#'     B-cell deconvolution. Only these genes are included in the fitting
#'     procedure to improve the robustness and biological specificity of the model.
#'   }
#' }
#'
#' @format
#' A list of length 3:
#' \itemize{
#'   \item \code{refProfiles}: double matrix, approximately \code{49902 × 6}
#'   \item \code{refProfiles.var}: double matrix, approximately \code{49902 × 6}
#'   \item \code{sigGenes}: character vector, approximately \code{65} genes
#' }
#'
#' @usage data(BRef)
#' @keywords datasets
"BRef"







#' Reference profiles for immune cell types using lm22 (EPIC/IOBR)
#'
#' \code{lm22} is a gene expression signature matrix used in the EPIC and IOBR
#' packages for estimating the proportions of different immune cell types from
#' bulk RNA-seq data. It contains reference profiles for multiple immune cell types,
#' including various B cell and T cell subsets, as well as plasma cells and
#' other immune cell types.
#'
#' @format
#' A list of length 2:
#' \describe{
#'   \item{\code{refProfiles}}{
#'     A numeric matrix where rows represent genes and columns represent different 
#'     immune cell types. The matrix provides gene expression profiles used as references 
#'     for deconvolution in the EPIC model. The number of columns may vary depending on the 
#'     number of immune cell types included in the reference profiles.
#'   }
#'   \item{\code{cellTypes}}{
#'     A character vector containing the names of the immune cell types included 
#'     in the reference profiles.
#'   }
#' }
#'
#' @usage data(lm22)
#' @keywords datasets
"lm22"



#' GSE62254 Gastric Cancer Expression Data
#'
#' Gene expression data from gastric cancer patients in the GSE62254 dataset.
#' Contains raw count matrix from RNA-seq or microarray experiments.
#'
#' @format A numeric matrix with genes as rows and samples as columns:
#' \describe{
#'   \item{Rows}{Ensembl gene identifiers (e.g., ENSG00000000003)}
#'   \item{Columns}{Sample identifiers (e.g., TCGA-2F-A9KQ)}
#'   \item{Values}{Expression counts or intensities}
#' }
#'
#' @usage data(eset_gse62254)
#'
#' @source GEO accession GSE62254
#'
#' @references
#' Cristescu R et al. (2015) Molecular analysis of gastric cancer identifies
#' subtypes associated with distinct clinical outcomes. Nat Med 21:449-456.
#' doi:10.1038/nm.3850
#'
#' @keywords datasets
"eset_gse62254"


#' TCGA-STAD Tumor Microenvironment Signature Scores
#'
#' Pre-calculated tumor microenvironment (TME) signature scores for TCGA
#' stomach adenocarcinoma (STAD) samples. Contains expression levels of
#' key TME-related genes and potentially immune/stomal signature scores.
#'
#' @format A numeric matrix with genes/signatures as rows and samples as columns:
#' \describe{
#'   \item{Rows}{Gene symbols or signature names (e.g., B2M, HLA-B, HSPB1)}
#'   \item{Columns}{TCGA sample identifiers (e.g., TCGA-3M-AB46-01A)}
#'   \item{Values}{Normalized expression values (appears to be log2-scale)}
#' }
#'
#' @usage data(eset_tme_stad)
#'
#' @source The Cancer Genome Atlas Stomach Adenocarcinoma (TCGA-STAD)
#'
#' @references
#' Cancer Genome Atlas Research Network. Comprehensive molecular
#' characterization of gastric adenocarcinoma. Nature. 2014;513:202-209.
#' doi:10.1038/nature13480
#'
#' @keywords datasets
"eset_tme_stad"


#' TCGA-BLCA Bladder Cancer Expression Data
#'
#' Gene expression data from bladder cancer patients in TCGA-BLCA.
#' Contains raw count matrix suitable for differential expression analysis.
#'
#' @format A numeric matrix with genes as rows and samples as columns:
#' \describe{
#'   \item{Rows}{Ensembl gene identifiers (e.g., ENSG00000000003)}
#'   \item{Columns}{TCGA sample identifiers (e.g., TCGA-2F-A9KO)}
#'   \item{Values}{Raw expression counts (integer values)}
#' }
#'
#' @usage data(eset_blca)
#'
#' @source The Cancer Genome Atlas Bladder Urothelial Carcinoma (TCGA-BLCA)
#'
#' @references
#' Cancer Genome Atlas Research Network. Comprehensive molecular
#' characterization of urothelial bladder carcinoma. Nature 507, 315-322 (2014).
#' doi:10.1038/nature12965
#'
#' @keywords datasets
"eset_blca"

#' TCGA-STAD Clinical and Molecular Annotation Data
#'
#' Clinical, molecular, and signature score data for TCGA stomach adenocarcinoma
#' (STAD) samples. Includes patient demographics, tumor characteristics,
#' molecular subtypes, and pre-computed signature scores.
#'
#'@format A data frame with samples as rows and variables as columns:
#' \itemize{
#'   \item ID – TCGA sample barcode
#'   \item stage – tumor stage (Stage_I, Stage_II, etc.)
#'   \item status – vital status (Alive/Dead)
#'   \item Lauren – histological classification
#'   \item subtype – molecular subtype (EBV, MSI, GS, CN, ...)
#'   \item EBV – EBV infection status
#'   \item TMEscore_plus – continuous tumor-micro-environment score
#'   \item TMEscore_plus_binary – High/Low TME classification
#'   \item time – follow-up time (months)
#'   \item OS_status – 0 = alive, 1 = dead
#'   \item ARID1A, PIK3CA – driver-gene mutation status
#'   \item MALAT1 – lncRNA expression
#'   \item remaining columns – gene-expression values and additional clinical/molecular annotations
#' }
#'
#'
#' @usage data(tcga_stad_pdata)
#'
#' @source The Cancer Genome Atlas Stomach Adenocarcinoma (TCGA-STAD)
#'
#' @references
#' Cancer Genome Atlas Research Network. Comprehensive molecular
#' characterization of gastric adenocarcinoma. Nature 513, 202-209 (2014).
#' doi:10.1038/nature13480
#'
#' @keywords datasets
"tcga_stad_pdata"


#' IMvigor210 Bladder Cancer Immunotherapy Cohort Data
#'
#' Clinical and biomarker data from the IMvigor210 clinical trial cohort.
#' Includes treatment response, survival outcomes, and immune biomarker
#' measurements for bladder cancer patients treated with atezolizumab.
#'
#' @format A data frame with patients as rows and variables as columns:
#' \describe{
#'   \item{ID}{Patient sample identifier}
#'   \item{BOR}{Best overall response (CR, PR, SD, PD, NA)}
#'   \item{BOR_binary}{Binary response classification (R=responder, NR=non-responder)}
#'   \item{OS_days}{Overall survival time in days}
#'   \item{OS_status}{Overall survival status (0=alive, 1=dead)}
#'   \item{Mutation_Load}{Tumor mutation burden}
#'   \item{Neo_antigen_Load}{Neoantigen load}
#'   \item{CD_8_T_effector}{CD8+ T effector signature score}
#'   \item{Immune_Checkpoint}{Immune checkpoint signature score}
#'   \item{Pan_F_TBRs}{Pan-fibroblast TGF-β response signature}
#'   \item{Mismatch_Repair}{Mismatch repair status or signature}
#'   \item{TumorPurity}{Estimated tumor purity}
#' }
#'
#' @usage data(imvigor210_pdata)
#'
#' @source IMvigor210 clinical trial (NCT02108652)
#'
#' @references
#' Mariathasan S et al. TGFβ attenuates tumour response to PD-L1 blockade
#' by contributing to exclusion of T cells. Nature 554, 544-548 (2018).
#' doi:10.1038/nature25501
#'
#' @keywords datasets
"imvigor210_pdata"


#' TCGA-STAD Gastric Cancer Cohort with Molecular and Clinical Data
#'
#' A comprehensive dataset containing clinical, pathological, and molecular data
#' from The Cancer Genome Atlas (TCGA) Stomach Adenocarcinoma (STAD) project.
#' Includes RNA-seq expression data, survival outcomes, and pathological features
#' for gastric cancer patients.
#'
#' @format A data frame with 374 rows (patients) and 323 variables:
#' \itemize{
#'   \item ID – TCGA barcode
#'   \item ProjectID – "TCGA-STAD"
#'   \item Technology / platform – sequencing details
#'   \item Gender – M/F
#'   \item Age – age at diagnosis (years)
#'   \item Survival outcomes – RFS_status, OS_time, OS_status
#'   \item Pathology – Lauren type, differentiation, AJCC_stage, T/N/M_stage
#'   \item 308 omics columns – gene-level RNA-seq counts / TPM / signatures
#' }
#'
#' 
#' @usage data(sig_stad)
#'
#' @source The Cancer Genome Atlas (TCGA) Research Network
#' \url{https://www.cancer.gov/tcga}
#'
#' @references
#' Cancer Genome Atlas Research Network. Comprehensive molecular characterization
#' of gastric adenocarcinoma. Nature 513, 202-209 (2014).
#' doi:10.1038/nature13480
#'
#' Liu J et al. Integrated omics analysis of gastric cancer.
#' Cell Reports 29, 1-15 (2019).
#' doi:10.1016/j.celrep.2019.09.045
#'
#' @keywords datasets
#' @docType data
"sig_stad"

#' TCGA-STAD Gastric Cancer Immune Infiltration Signatures
#'
#' A comprehensive dataset containing immune cell infiltration scores derived
#' from CIBERSORT analysis of TCGA Stomach Adenocarcinoma (STAD) RNA-seq data.
#' Includes detailed immune cell subtype proportions for each gastric cancer sample.
#'
#'@format A data frame with 375 rows (samples) and 464 variables:
#' \itemize{
#'   \item ID – TCGA barcode
#'   \item 22 CIBERSORT immune-cell proportions (0–1 scale)
#'   \item Immune-checkpoint and signature scores
#'   \item Remaining columns – other deconvolution / signature outputs
#' }
#' 
#' @usage data(tcga_stad_sig)
#'
#' @source The Cancer Genome Atlas (TCGA) Research Network, processed with CIBERSORT
#' \url{https://www.cancer.gov/tcga}
#' \url{https://cibersort.stanford.edu/}
#'
#' @references
#' Newman AM et al. Robust enumeration of cell subsets from tissue expression profiles.
#' Nature Methods 12, 453-457 (2015).
#' doi:10.1038/nmeth.3337
#'
#' Thorsson V et al. The Immune Landscape of Cancer.
#' Immunity 48, 812-830 (2018).
#' doi:10.1016/j.immuni.2018.03.023
#'
#' Cancer Genome Atlas Research Network. Comprehensive molecular characterization
#' of gastric adenocarcinoma. Nature 513, 202-209 (2014).
#' doi:10.1038/nature13480
#'
#' @keywords datasets
#' @docType data
"tcga_stad_sig"


#' IMvigor210 Bladder Cancer Cohort Multi-omics Signatures
#'
#' Comprehensive multi-omics dataset from the IMvigor210 phase II clinical trial
#' of metastatic urothelial cancer patients treated with atezolizumab (anti-PD-L1).
#' This dataset includes CIBERSORT immune cell deconvolution scores, gene expression
#' values, and various immune-related molecular signatures.
#'
#' @format A data frame with 348 rows (patients) and 456 variables:
#' \itemize{
#'   \item ID – unique sample identifier ("SAM***")
#'   \item 22 CIBERSORT immune-cell proportions (0–1 scale)
#'   \item immune-gene expression values
#'   \item immune signature scores (GEP, cytolytic activity, etc.)
#'   \item clinical / molecular features extracted from the original publication
#' }
#' 
#' @usage data(imvigor210_sig)
#'
#' @source IMvigor210 clinical trial (NCT02108652)
#' Supplementary data from: Mariathasan S et al. Nature 554, 544-548 (2018)
#' \url{https://www.nature.com/articles/nature25501}
#' \url{https://clinicaltrials.gov/ct2/show/NCT02108652}
#'
#' @references
#' Mariathasan S et al. TGFβ attenuates tumour response to PD-L1 blockade
#' by contributing to exclusion of T cells. Nature 554, 544-548 (2018).
#' doi:10.1038/nature25501
#'
#' Rosenberg JE et al. Atezolizumab in patients with locally advanced and
#' metastatic urothelial carcinoma who have progressed following treatment with
#' platinum-based chemotherapy: a single-arm, multicentre, phase 2 trial.
#' Lancet 387, 1909-1920 (2016). doi:10.1016/S0140-6736(16)00561-4
#'
#' Newman AM et al. Robust enumeration of cell subsets from tissue expression profiles.
#' Nature Methods 12, 453-457 (2015). doi:10.1038/nmeth.3337
#'
#' @keywords datasets
#' @docType data
"imvigor210_sig"


#' Single-cell RNA-seq Differential Expression Analysis Results
#'
#' Example dataset containing differential expression analysis results
#' from single-cell RNA sequencing (scRNA-seq) analysis.
#' This dataset serves as input for signature analysis functions in the IOBR package,
#' particularly \code{\link{get_sig_sc}} and \code{\link{sig_gsea}}.
#'
#' @format A data frame with 3,212 rows (genes) and 7 columns (statistics):
#' \describe{
#'   \item{p_val}{Raw p-value for differential expression test}
#'   \item{avg_log2FC}{Average log2 fold-change of gene expression}
#'   \item{pct.1}{Percentage of cells expressing the gene in cluster 1}
#'   \item{pct.2}{Percentage of cells expressing the gene in cluster 2}
#'   \item{p_val_adj}{Adjusted p-value (e.g., Bonferroni, BH FDR correction)}
#'   \item{cluster}{Cell cluster or cell type identifier (e.g., "Epithelial cells 2")}
#'   \item{gene}{Gene symbol (e.g., "IGFBP3", "PCDH7")}
#' }
#'
#' @usage data(deg)
#'
#' @keywords datasets
#' @docType data
"deg"

#' Example Dataset for Subgroup Survival Analysis
#'
#' An example dataset demonstrating the data structure required for subgroup
#' survival analysis using the \code{\link{subgroup_survival}} function.
#' Contains simulated clinical and biomarker data with survival outcomes.
#'
#' @format A data frame with clinical variables and biomarker scores:
#' \describe{
#'   \item{Patient_ID}{Unique patient identifier}
#'   \item{ProjectID}{Study or dataset identifier (e.g., "Dataset1")}
#'   \item{AJCC_stage}{AJCC pathological stage (2, 3, 4)}
#'   \item{status}{Event status (0=censored, 1=event)}
#'   \item{time}{Follow-up time in months}
#'   \item{score}{Continuous biomarker score (numeric values)}
#'   \item{score_binary}{Binary biomarker classification ("High", "Low")}
#' }
#'
#' @usage data(subgroup_data)
#'
#' @keywords datasets
#' @docType data
"subgroup_data"

#' Example Clinical Data for TCGA-STAD Gastric Cancer Analysis
#'
#' A small example dataset demonstrating the clinical data structure required
#' for TCGA Stomach Adenocarcinoma (STAD) analysis using IOBR package functions.
#' Contains simulated clinical variables, molecular subtypes, and survival data.
#'
#' @format A data frame with patient samples as rows and clinical variables as columns:
#' \describe{
#'   \item{ID}{Unique patient identifier (TCGA barcode format)}
#'   \item{stage}{AJCC pathological stage (Stage_II, Stage_III, Stage_IV)}
#'   \item{status}{Patient vital status (Alive, Dead, NA)}
#'   \item{Lauren}{Lauren histological classification (Intestinal, Diffuse, Mixed)}
#'   \item{subtype}{Molecular subtype classification (EBV, GS)}
#'   \item{EBV}{Epstein-Barr virus status (Positive, Negative)}
#'   \item{time}{Overall survival time in months}
#'   \item{OS_status}{Overall survival status (0=alive, 1=dead)}
#' }
#'
#' @usage data(stad_group)
#'
#' @keywords datasets
#' @docType data
"stad_group"




#' GRCh38 Human Genome Annotation
#'
#' Gene annotation for human GRCh38/hg38 genome assembly, used for expression
#' normalization and gene identifier mapping in IOBR functions.
#'
#' @format A data frame with 60668 rows and 11 columns:
#' \describe{
#'   \item{id}{Ensembl gene identifier.}
#'   \item{eff_length}{Effective gene length for TPM calculation.}
#'   \item{gc}{GC content proportion.}
#'   \item{entrez}{Entrez Gene ID.}
#'   \item{symbol}{Official gene symbol.}
#'   \item{chr}{Chromosome name.}
#'   \item{start}{Genomic start position.}
#'   \item{end}{Genomic end position.}
#'   \item{strand}{Genomic strand.}
#'   \item{biotype}{Gene biotype classification.}
#'   \item{description}{Gene description.}
#' }
#'
#' @usage data(anno_grch38)
#'
#' @source Ensembl release 104 (GRCh38.p13)
#'
#' @keywords datasets
"anno_grch38"


#' Mouse Genome Annotation (GC/VM32)
#'
#' Mouse gene annotation dataset for genome assembly GRCm38/mm10,
#' containing gene features including GC content and gene length information.
#' Used for mouse transcriptomic data analysis.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{id}{Ensembl gene identifier (e.g., "ENSMUSG00000000001").}
#'   \item{eff_length}{Effective gene length for TPM calculation.}
#'   \item{gc}{GC content proportion (0-1 range).}
#'   \item{symbol}{Official gene symbol.}
#'   \item{mgi_id}{Mouse Genome Informatics identifier.}
#'   \item{gene_type}{Gene type classification (e.g., "protein_coding", "lncRNA").}
#'   \item{start}{Genomic start position.}
#'   \item{end}{Genomic end position.}
#'   \item{transcript_id}{Transcript identifier (mostly NA in this dataset).}
#'   \item{ont}{Gene ontology information (mostly NA in this dataset).}
#' }
#'
#' @usage data(anno_gc_vm32)
#'
#' @source Ensembl database for mouse GRCm38/mm10
#'
#' @keywords datasets
"anno_gc_vm32"


#' Affymetrix Human Genome U133 Plus 2.0 Array Annotation
#'
#' Probe annotation for the Affymetrix Human Genome U133 Plus 2.0 microarray.
#' Provides mapping between probe IDs and gene symbols for microarray data
#' analysis.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{probe_id}{Affymetrix probe identifier.}
#'   \item{symbol}{Official gene symbol.}
#' }
#'
#' @usage data(anno_hug133plus2)
#'
#' @source Affymetrix annotation files (HuGene-1_0-st-v1)
#'
#' @keywords datasets
"anno_hug133plus2"

#' General RNA-seq Annotation
#'
#' Generic gene annotation for RNA-seq data analysis. Provides basic
#' gene identifier mappings for various RNA-seq analysis workflows.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{probe_id}{Gene identifier (platform specific).}
#'   \item{symbol}{Official gene symbol.}
#' }
#'
#' @usage data(anno_rnaseq)
#'
#' @source Curated from multiple public RNA-seq resources
#'
#' @keywords datasets
"anno_rnaseq"

#' Illumina Microarray Annotation
#'
#' Probe annotation for Illumina microarray platforms. Provides mapping
#' between Illumina probe IDs and gene symbols for microarray data analysis.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{probe_id}{Illumina probe identifier.}
#'   \item{symbol}{Official gene symbol.}
#' }
#'
#' @usage data(anno_illumina)
#'
#' @source Illumina annotation manifest files
#'
#' @keywords datasets
"anno_illumina"



