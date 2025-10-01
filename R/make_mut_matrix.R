#' Construct Mutation Matrices from MAF Data
#'
#' Builds mutation presence/absence matrices from MAF input (file path or MAF object). Supports multiple categories: all mutations, SNPs, indels, and frameshift mutations. When category = "multi", returns a list of matrices for each category. Compatible with TCGA-formatted data.
#'
#' @param maf Character or MAF object. Path to MAF file or an already loaded MAF object.
#' @param mut_data Data frame or NULL. Preloaded MAF-like data (used if 'maf' is NULL).
#' @param isTCGA Logical. Whether the MAF follows TCGA conventions. Default is TRUE.
#' @param category Character. Mutation category: "all", "snp", "indel", "frameshift", or "multi". Default is "multi".
#' @param Tumor_Sample_Barcode Character. Column name for tumor sample IDs. Default is "Tumor_Sample_Barcode".
#' @param Hugo_Symbol Character. Column name for gene symbols. Default is "Hugo_Symbol".
#' @param Variant_Classification Character. Column name for variant classification (e.g., Frame_Shift_Del). Default is "Variant_Classification".
#' @param Variant_Type Character. Column name for variant type (e.g., SNP, INS, DEL). Default is "Variant_Type".
#'
#' @return List of mutation matrices (if category = "multi") or a single matrix for the specified category.
#' @export
#' @author Dongqian Zeng
#' @author Shixiang Huang
#'
#' @examples
#' # See maftools and TCGAbiolinks documentation for obtaining MAF input
#' # mut_list <- make_mut_matrix(maf = maf, isTCGA = TRUE, category = "multi")
make_mut_matrix <- function(maf = NULL, mut_data = NULL, isTCGA = TRUE, category = "multi", Tumor_Sample_Barcode = "Tumor_Sample_Barcode", Hugo_Symbol = "Hugo_Symbol", Variant_Classification = "Variant_Classification", Variant_Type = "Variant_Type") {
  if (!is.null(maf)) {
    mut_maf <- maftools::read.maf(maf = maf, useAll = TRUE, isTCGA = isTCGA)

    print(summary(mut_maf@data$Variant_Classification))
    print(summary(mut_maf@data$Variant_Type))
    #########################################

    mut <- mut_maf@data
  } else {
    mut <- as.data.frame(mut_data)
    colnames(mut)[which(colnames(mut) == Tumor_Sample_Barcode)] <- "Tumor_Sample_Barcode"
    colnames(mut)[which(colnames(mut) == Hugo_Symbol)] <- "Hugo_Symbol"
    colnames(mut)[which(colnames(mut) == Variant_Classification)] <- "Variant_Classification"

    if (!Variant_Type %in% colnames(mut)) {
      if ("filter" %in% colnames(mut)) mut <- mut[mut$filter == "PASS", ]
      mut$Variant_Type <- mut$Variant_Classification
      mut <- mut[!grepl(mut$Variant_Type, pattern = "synonymous"), ]

      mut$Variant_Type[stringr::str_detect(tolower(mut$Variant_Type), pattern = "missense")] <- "SNP"
      mut$Variant_Type[stringr::str_detect(tolower(mut$Variant_Type), pattern = "frameshift")] <- "Frame_Shift"
      mut$Variant_Type[stringr::str_detect(tolower(mut$Variant_Type), pattern = "insert")] <- "INS"
      mut$Variant_Type[stringr::str_detect(tolower(mut$Variant_Type), pattern = "delet")] <- "DEL"

      mut <- mut[mut$Variant_Type %in% c("SNP", "INS", "DEL", "Frame_Shift"), ]
    } else {
      colnames(mut)[which(colnames(mut) == Variant_Type)] <- "Variant_Type"
    }
  }

  mut <- mut[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Variant_Type")]

  if (category == "multi") {
    mut_all <- reshape2::dcast(mut, Hugo_Symbol ~ Tumor_Sample_Barcode, value.var = "Variant_Classification")
    mut_all <- tibble::column_to_rownames(mut_all, var = "Hugo_Symbol")
    mut_all <- as.matrix(mut_all)
    mut_all <- matrix(as.numeric(mut_all), dim(mut_all), dimnames = dimnames(mut_all))

    mut_snp <- mut[mut$Variant_Type == "SNP", ]
    mut_snp <- reshape2::dcast(mut_snp, Hugo_Symbol ~ Tumor_Sample_Barcode, value.var = "Variant_Classification")
    mut_snp <- tibble::column_to_rownames(mut_snp, var = "Hugo_Symbol")
    mut_snp <- as.matrix(mut_snp)
    mut_snp <- matrix(as.numeric(mut_snp), dim(mut_snp), dimnames = dimnames(mut_snp))

    mut_indel <- mut[mut$Variant_Type == "DEL" | mut$Variant_Type == "INS", ]
    mut_indel <- reshape2::dcast(mut_indel, Hugo_Symbol ~ Tumor_Sample_Barcode, value.var = "Variant_Classification")
    mut_indel <- tibble::column_to_rownames(mut_indel, var = "Hugo_Symbol")
    mut_indel <- as.matrix(mut_indel)
    mut_indel <- matrix(as.numeric(mut_indel), dim(mut_indel), dimnames = dimnames(mut_indel))

    mut_frameshift <- mut[grepl(tolower(mut$Variant_Classification), pattern = "frame"), ]
    mut_frameshift <- reshape2::dcast(mut_frameshift, Hugo_Symbol ~ Tumor_Sample_Barcode, value.var = "Variant_Classification")
    mut_frameshift <- tibble::column_to_rownames(mut_frameshift, var = "Hugo_Symbol")
    mut_frameshift <- as.matrix(mut_frameshift)
    mut_frameshift <- matrix(as.numeric(mut_frameshift), dim(mut_frameshift), dimnames = dimnames(mut_frameshift))
    ####################################
    ####################################
    mut_list <- list(all = t(mut_all), snp = t(mut_snp), indel = t(mut_indel), frameshift = t(mut_frameshift))
    mut_list <- ComplexHeatmap::unify_mat_list(mut_list)
  } else {
    mut <- mut[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Variant_Type")]
    mut_all <- reshape2::dcast(mut, Hugo_Symbol ~ Tumor_Sample_Barcode, value.var = "Variant_Classification")
    mut_all <- tibble::column_to_rownames(mut_all, var = "Hugo_Symbol")
    mut_all <- as.matrix(mut_all)
    mut_all <- matrix(as.numeric(mut_all), dim(mut_all), dimnames = dimnames(mut_all))
    mut_list <- list(all = t(mut_all))
  }
  return(mut_list)
}
