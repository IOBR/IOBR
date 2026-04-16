#' Construct Mutation Matrices from MAF Data
#'
#' Builds mutation presence/absence matrices from MAF input (file path or MAF
#' object). Supports multiple categories: all mutations, SNPs, indels, and
#' frameshift mutations. When category = "multi", returns a list of matrices
#' for each category. Compatible with TCGA-formatted data.
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
#' @note
#' Some users may encounter errors from upstream data import (e.g.
#' "Can't combine ..$Tumor_Seq_Allele2" when using TCGAbiolinks or
#' TCGAmutations).
#' This is due to inconsistent column types in the source MAF tables,
#' not an issue of this function.
#' Please ensure your MAF or merged data frame uses consistent column types
#' (e.g. convert allele columns to character before input).
#'
#' @export
#' @author Dongqian Zeng
#' @author Shixiang Huang
#'
#' @examples
#' \dontrun{
#' # See maftools or TCGAbiolinks documentation for obtaining MAF input
#' # Example: Download MAF file from TCGA portal
#' mut_list <- make_mut_matrix(maf = "path_to_maf_file.maf", isTCGA = TRUE, category = "multi")
#' }
make_mut_matrix <- function(maf = NULL, mut_data = NULL, isTCGA = TRUE,
                            category = c("multi", "all", "snp", "indel", "frameshift"),
                            Tumor_Sample_Barcode = "Tumor_Sample_Barcode",
                            Hugo_Symbol = "Hugo_Symbol",
                            Variant_Classification = "Variant_Classification",
                            Variant_Type = "Variant_Type") {
  category <- rlang::arg_match(category)

  if (is.null(maf) && is.null(mut_data)) {
    cli::cli_abort("Either {.arg maf} or {.arg mut_data} must be provided.")
  }

  if (!is.null(maf)) {
    mut_maf <- maftools::read.maf(maf = maf, useAll = TRUE, isTCGA = isTCGA)

    cli::cli_alert_info("Variant Classification summary:")
    print(summary(mut_maf@data$Variant_Classification))
    cli::cli_alert_info("Variant Type summary:")
    print(summary(mut_maf@data$Variant_Type))

    mut <- mut_maf@data
  } else {
    mut <- as.data.frame(mut_data)

    required_cols <- c(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification)
    missing_cols <- setdiff(required_cols, colnames(mut))
    if (length(missing_cols) > 0) {
      cli::cli_abort("Missing required column(s): {.val {missing_cols}}")
    }

    colnames(mut)[colnames(mut) == Tumor_Sample_Barcode] <- "Tumor_Sample_Barcode"
    colnames(mut)[colnames(mut) == Hugo_Symbol] <- "Hugo_Symbol"
    colnames(mut)[colnames(mut) == Variant_Classification] <- "Variant_Classification"

    if (!Variant_Type %in% colnames(mut)) {
      if ("filter" %in% colnames(mut)) mut <- mut[mut$filter == "PASS", ]
      mut$Variant_Type <- mut$Variant_Classification
      mut <- mut[!grepl("synonymous", mut$Variant_Type, ignore.case = TRUE), ]

      mut$Variant_Type <- dplyr::case_when(
        stringr::str_detect(tolower(mut$Variant_Type), "missense") ~ "SNP",
        stringr::str_detect(tolower(mut$Variant_Type), "frameshift") ~ "Frame_Shift",
        stringr::str_detect(tolower(mut$Variant_Type), "insert") ~ "INS",
        stringr::str_detect(tolower(mut$Variant_Type), "delet") ~ "DEL",
        TRUE ~ mut$Variant_Type
      )

      mut <- mut[mut$Variant_Type %in% c("SNP", "INS", "DEL", "Frame_Shift"), ]
    } else {
      colnames(mut)[colnames(mut) == Variant_Type] <- "Variant_Type"
    }
  }

  mut <- mut[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Variant_Type")]

  mut_list <- .build_mut_matrices(mut, category)

  return(mut_list)
}

#' Build Mutation Matrices from MAF Data
#'
#' @param mut Data frame with mutation data.
#' @param category Character. Mutation category.
#'
#' @return List of mutation matrices.
#'
#' @keywords internal
.build_mut_matrices <- function(mut, category) {
  if (category == "multi") {
    mut_all <- .create_mut_matrix(mut)
    mut_snp <- .create_mut_matrix(mut[mut$Variant_Type == "SNP", ])
    mut_indel <- .create_mut_matrix(mut[mut$Variant_Type %in% c("DEL", "INS"), ])
    mut_frameshift <- .create_mut_matrix(mut[grepl("frame", tolower(mut$Variant_Classification)), ])

    mut_list <- list(
      all = t(mut_all),
      snp = t(mut_snp),
      indel = t(mut_indel),
      frameshift = t(mut_frameshift)
    )
    rlang::check_installed("ComplexHeatmap")
    mut_list <- ComplexHeatmap::unify_mat_list(mut_list)
  } else {
    mut_all <- .create_mut_matrix(mut)
    mut_list <- list(all = t(mut_all))
  }

  mut_list
}

#' Create Single Mutation Matrix
#'
#' @param mut Data frame with mutation data.
#'
#' @return Mutation matrix.
#'
#' @keywords internal
.create_mut_matrix <- function(mut) {
  if (nrow(mut) == 0) {
    cli::cli_warn("No mutations found for the specified category.")
    return(matrix(nrow = 0, ncol = 0))
  }

  rlang::check_installed("reshape2")
  mut_mat <- reshape2::dcast(mut, Hugo_Symbol ~ Tumor_Sample_Barcode,
    value.var = "Variant_Classification"
  )
  mut_mat <- tibble::column_to_rownames(mut_mat, var = "Hugo_Symbol")
  mut_mat <- as.matrix(mut_mat)
  mut_mat <- matrix(as.numeric(mut_mat),
    nrow = nrow(mut_mat), ncol = ncol(mut_mat),
    dimnames = dimnames(mut_mat)
  )
  mut_mat
}
