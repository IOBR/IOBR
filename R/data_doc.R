#' NULL Model Coefficients for MCPcounter
#'
#' @format A `data.frame` with cell types in rows and coefficients in columns.
#' @usage data(null_models)
#' @keywords datasets
#'
#' @examples
#' data(null_models)
#' head(null_models)
"null_models"


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
#'   \item{time}{Numeric. Overall survival or follow-up time, typically
#'   measured in months.}
#'   \item{OS_status}{Integer/binary. Overall survival status indicator.
#'   (\code{1} = death, \code{0} = censored)}
#' }
#'
#' @usage data(pdata_stad)
#' @keywords datasets
#'
#' @examples
#' data(pdata_stad)
#' head(pdata_stad)
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
#'
#' @examples
#' data(signature_collection)
#' head(signature_collection)
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
#'     \item{EMT}{Epithelial–mesenchymal transition (EMT)–associated
#'     signatures.}
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
#'
#' @examples
#' data(sig_group)
#' head(sig_group)
"sig_group"

#' TCGA-STAD Clinical and Molecular Annotation Data
#'
#' Clinical, molecular, and signature score data for TCGA stomach adenocarcinoma
#' (STAD) samples. Includes patient demographics, tumor characteristics,
#' molecular subtypes, and pre-computed signature scores.
#'
#' @format A data frame with samples as rows and variables as columns:
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
#'   \item remaining columns – gene-expression values and additional
#'   clinical/molecular annotations
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
#'
#' @examples
#' data(tcga_stad_pdata)
#' head(tcga_stad_pdata)
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
#'   \item{BOR_binary}{Binary response classification (R=responder,
#'   NR=non-responder)}
#'   \item{OS_days}{Overall survival time in days}
#'   \item{OS_status}{Overall survival status (0=alive, 1=dead)}
#'   \item{Mutation_Load}{Tumor mutation burden}
#'   \item{Neo_antigen_Load}{Neoantigen load}
#'   \item{CD_8_T_effector}{CD8+ T effector signature score}
#'   \item{Immune_Checkpoint}{Immune checkpoint signature score}
#'   \item{Pan_F_TBRs}{Pan-fibroblast TGF-\eqn{\beta} response signature}
#'   \item{Mismatch_Repair}{Mismatch repair status or signature}
#'   \item{TumorPurity}{Estimated tumor purity}
#' }
#'
#' @usage data(imvigor210_pdata)
#'
#' @source IMvigor210 clinical trial (NCT02108652)
#'
#' @references
#' Mariathasan S et al. TGF\eqn{\beta} attenuates tumour response to PD-L1 blockade
#' by contributing to exclusion of T cells. Nature 554, 544-548 (2018).
#' doi:10.1038/nature25501
#'
#' @keywords datasets
#'
#' @examples
#' data(imvigor210_pdata)
#' head(imvigor210_pdata)
"imvigor210_pdata"




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
#'
#' @examples
#' data(subgroup_data)
#' head(subgroup_data)
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
#'   \item{Lauren}{Lauren histological classification (Intestinal, Diffuse,
#'   Mixed)}
#'   \item{subtype}{Molecular subtype classification (EBV, GS)}
#'   \item{EBV}{Epstein-Barr virus status (Positive, Negative)}
#'   \item{time}{Overall survival time in months}
#'   \item{OS_status}{Overall survival status (0=alive, 1=dead)}
#' }
#'
#' @usage data(stad_group)
#'
#' @keywords datasets
#'
#' @examples
#' data(stad_group)
#' head(stad_group)
#' @docType data
"stad_group"




