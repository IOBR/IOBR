#' Convert Read Counts to Transcripts Per Million (TPM)
#'
#' @description
#' Transforms gene expression count data into Transcripts Per Million (TPM)
#' values, normalizing for gene length and library size. Supports multiple gene
#' ID types and can retrieve gene length information from BioMart or use local
#' datasets.
#'
#' @param countMat Numeric matrix of raw read counts with genes in rows and
#'   samples in columns.
#' @param idType Character string specifying the gene identifier type. Options
#'   are `"Ensembl"`, `"Entrez"`, or `"Symbol"`. Default is `"Ensembl"`.
#' @param org Character string specifying the organism. Options include
#'   `"hsa"` (human) or `"mmus"` (mouse). Default is `"hsa"`.
#' @param source Character string specifying the source for gene length
#'   information. Options are `"biomart"` (retrieve from Ensembl BioMart) or
#'   `"local"` (use local dataset). Default is `"local"`.
#' @param effLength Data frame containing effective gene length information. If
#'   `NULL`, lengths are retrieved based on `source`. Default is `NULL`.
#' @param id Character string specifying the column name in `effLength`
#'   containing gene identifiers. Default is `"id"`.
#' @param gene_symbol Character string specifying the column name in `effLength`
#'   containing gene symbols. Default is `"symbol"`.
#' @param length Character string specifying the column name in `effLength`
#'   containing gene lengths. Default is `"eff_length"`.
#' @param check_data Logical indicating whether to check for missing values in
#'   the count matrix. Default is `FALSE`.
#'
#' @return Data frame of TPM-normalized expression values with genes in rows
#'   and samples in columns.
#'
#' @author Wubing Zhang, Dongqiang Zeng, Yiran Fang
#' @export
#'
#' @examples
#' \donttest{
#' # Load TCGA count data
#' eset_stad <- load_data("eset_stad")
#'
#' # Transform to TPM using local gene annotation
#' eset <- count2tpm(countMat = eset_stad, source = "local", idType = "ensembl")
#' head(eset)
#' }
count2tpm <- function(countMat,
                      idType = "Ensembl",
                      org = c("hsa", "mmus"),
                      source = c("local", "biomart"),
                      effLength = NULL,
                      id = "id",
                      gene_symbol = "symbol",
                      length = "eff_length",
                      check_data = FALSE) {
  # Validate arguments
  org <- rlang::arg_match(org)
  source <- rlang::arg_match(source)

  # Ensure matrix format
  if (!is.matrix(countMat)) {
    countMat <- as.matrix(countMat)
  }
  storage.mode(countMat) <- "numeric"

  # Check for missing values
  if (check_data || sum(is.na(countMat)) > 0) {
    na_count <- sum(is.na(countMat))
    if (na_count > 0) {
      cli::cli_alert_warning(
        "Found {na_count} missing value{?s}; removing affected genes"
      )
      feas <- feature_manipulation(data = countMat, is_matrix = TRUE)
      countMat <- countMat[rownames(countMat) %in% feas, , drop = FALSE]
    }
  }

  # Get gene lengths
  if (!is.null(effLength)) {
    # Use provided effective lengths
    len <- .get_lengths_from_df(countMat, effLength, id, gene_symbol, length)
  } else if (source == "biomart") {
    # Fetch from BioMart
    len <- .get_lengths_from_biomart(countMat, idType, org)
  } else {
    # Use local annotation
    len <- .get_lengths_local(countMat, idType, org)
  }

  # Remove genes without length information
  valid_idx <- !is.na(len)
  if (!all(valid_idx)) {
    n_omit <- sum(!valid_idx)
    cli::cli_alert_warning("Omitting {n_omit} genes without length information")
    countMat <- countMat[valid_idx, , drop = FALSE]
    len <- len[valid_idx]
  }

  # Calculate TPM
  TPM <- .calculate_tpm(countMat, len)

  # Remove duplicates and clean up
  TPM <- TPM[!is.na(rownames(TPM)) & !rownames(TPM) == " ", , drop = FALSE]

  symbol.id <- rownames(TPM)
  TPM <- as.data.frame(TPM)
  TPM$symbol <- symbol.id

  TPM <- remove_duplicate_genes(eset = TPM, column_of_symbol = "symbol")

  as.data.frame(TPM)
}

# Helper: Calculate TPM values
.calculate_tpm <- function(counts, lengths) {
  # RPK: reads per kilobase
  rpk <- counts / (lengths / 1000)
  # TPM: transcripts per million
  tpm <- 1e6 * t(t(rpk) / colSums(rpk, na.rm = TRUE))
  tpm
}

# Helper: Get lengths from data frame
.get_lengths_from_df <- function(countMat, effLength, id, gene_symbol, length) {
  effLength <- as.data.frame(effLength)
  colnames(effLength)[colnames(effLength) == id] <- "id"
  colnames(effLength)[colnames(effLength) == length] <- "eff_length"
  effLength <- effLength[!duplicated(effLength$id), ]

  countMat <- countMat[rownames(countMat) %in% effLength$id, , drop = FALSE]

  if (nrow(countMat) == 0) {
    cli::cli_abort("No matching identifiers found in effLength data")
  }

  effLength <- effLength[effLength$id %in% rownames(countMat), ]

  if (id != gene_symbol) {
    colnames(effLength)[colnames(effLength) == gene_symbol] <- "gene_symbol"
  } else {
    effLength$gene_symbol <- effLength$id
  }

  rownames(countMat) <- effLength[
    match(rownames(countMat), effLength$id),
    "gene_symbol"
  ]

  effLength[match(rownames(countMat), effLength$gene_symbol), "eff_length"]
}

# Helper: Get lengths from BioMart
.get_lengths_from_biomart <- function(countMat, idType, org) {
  cli::cli_alert_warning(
    "BioMart source is being deprecated. Consider using source='local'"
  )

  rlang::check_installed("biomaRt")

  datasets <- paste0(c(
    "hsapiens", "mmusculus", "btaurus", "cfamiliaris",
    "ptroglodytes", "rnorvegicus", "sscrofa"
  ), "_gene_ensembl")

  type <- c(
    "ensembl_gene_id", "entrezgene_id", "hgnc_symbol",
    "start_position", "end_position"
  )
  if (org == "mmus") type[3] <- "mgi_symbol"

  ds <- datasets[grepl(org, datasets)]
  mart <- biomaRt::useMart(
    host = "https://www.ensembl.org",
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = ds
  )

  ensembl <- biomaRt::getBM(attributes = type, mart = mart)
  ensembl$Length <- abs(ensembl$end_position - ensembl$start_position)

  id_col <- switch(toupper(idType),
    "ENSEMBL" = "ensembl_gene_id",
    "ENTREZ"  = "entrezgene_id",
    "SYMBOL"  = type[3],
    cli::cli_abort("Invalid idType: {.val {idType}}")
  )

  # Update rownames and return lengths
  idx <- match(rownames(countMat), ensembl[[id_col]])
  rownames(countMat) <- ensembl[idx, type[3]]
  ensembl$Length[idx]
}

# Helper: Get lengths from local annotation
.get_lengths_local <- function(countMat, idType, org) {
  idType_lower <- tolower(idType)

  if (org == "hsa") {
    .get_lengths_human(countMat, idType_lower)
  } else if (org == "mmus") {
    .get_lengths_mouse(countMat, idType_lower)
  }
}

# Helper: Get lengths for human genes
.get_lengths_human <- function(countMat, idType) {
  anno_grch38 <- load_data("anno_grch38")

  cli::cli_alert_info(
    "Using local annotation (anno_grch38) for TPM conversion"
  )

  if (idType == "ensembl") {
    rownames(countMat) <- substring(rownames(countMat), 1, 15)
    length_df <- anno_grch38[, c("id", "eff_length", "symbol")]
  } else if (idType == "entrez") {
    cli::cli_alert_warning(
      "Entrez IDs may yield fuzzy results. Consider using Ensembl IDs"
    )
    length_df <- anno_grch38[, c("entrez", "eff_length", "symbol")]
    colnames(length_df)[1] <- "id"
    length_df <- length_df[!duplicated(length_df$id), ]
  } else if (idType == "symbol") {
    cli::cli_alert_warning(
      "Gene symbols may yield fuzzy results. Consider using Ensembl IDs"
    )
    length_df <- anno_grch38[, c("symbol", "eff_length", "gc")]
    colnames(length_df)[1] <- "id"
    length_df <- length_df[!duplicated(length_df$id), ]
  } else {
    cli::cli_abort("Invalid idType for local source: {.val {idType}}")
  }

  .match_and_extract_lengths(countMat, length_df)
}

# Helper: Get lengths for mouse genes
.get_lengths_mouse <- function(countMat, idType) {
  anno_gc_vm32 <- load_data("anno_gc_vm32")

  cli::cli_alert_info(
    "Using local annotation (anno_gc_vm32) for TPM conversion"
  )

  if (idType == "ensembl") {
    length_df <- anno_gc_vm32[, c("id", "eff_length", "symbol")]
  } else if (idType == "mgi") {
    cli::cli_alert_warning(
      "MGI IDs may yield fuzzy results. Consider using Ensembl IDs"
    )
    length_df <- anno_gc_vm32[, c("mgi_id", "eff_length", "symbol")]
    colnames(length_df)[1] <- "id"
    length_df <- length_df[!duplicated(length_df$id), ]
  } else if (idType == "symbol") {
    cli::cli_alert_warning(
      "Gene symbols may yield fuzzy results. Consider using Ensembl IDs"
    )
    length_df <- anno_gc_vm32[, c("symbol", "eff_length", "gc")]
    colnames(length_df)[1] <- "id"
    length_df <- length_df[!duplicated(length_df$id), ]
  } else {
    cli::cli_abort("Invalid idType for mouse: {.val {idType}}")
  }

  .match_and_extract_lengths(countMat, length_df)
}

# Helper: Match identifiers and extract lengths
.match_and_extract_lengths <- function(countMat, length_df) {
  # Sort by length (descending) for duplicate handling
  length_df <- length_df[order(length_df$eff_length, decreasing = TRUE), ]

  # Filter to common genes
  common_genes <- intersect(rownames(countMat), length_df$id)

  if (length(common_genes) == 0) {
    cli::cli_abort(
      paste("No matching identifiers found between count matrix and annotation")
    )
  }

  countMat <- countMat[rownames(countMat) %in% common_genes, , drop = FALSE]
  length_df <- length_df[length_df$id %in% common_genes, ]

  # Match and return lengths
  idx <- match(rownames(countMat), length_df$id)
  lengths <- length_df$eff_length[idx]

  # Update rownames to symbols
  rownames(countMat) <- length_df[idx, 3]

  lengths
}
