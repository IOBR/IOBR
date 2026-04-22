#' Transform Signature Data into List Format
#'
#' @description
#' Converts signature data from a data frame (with signatures as columns and
#' genes as rows) into a list format suitable for IOBR functions. Handles
#' NA values appropriately.
#'
#' @param sig_data Data frame with signature names as columns and genes in
#'   rows. Use `NA` for missing values.
#' @param save_signature Logical. Whether to save the signature list as RData.
#'   Default is `FALSE`.
#' @param output_path Character. Full path (without extension) for output file.
#'   Required if `save_signature = TRUE`. Default is `NULL`.
#'
#' @return List of signatures.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' sig_data <- data.frame(
#'   Signature1 = c("Gene1", "Gene2", "Gene3", NA),
#'   Signature2 = c("Gene4", "Gene5", NA, NA)
#' )
#' format_signatures(sig_data)
format_signatures <- function(sig_data, save_signature = FALSE, output_path = NULL) {
  if (!is.data.frame(sig_data) && !is.matrix(sig_data)) {
    cli::cli_abort("{.arg sig_data} must be a data frame or matrix")
  }
  if (ncol(sig_data) == 0) {
    cli::cli_abort("{.arg sig_data} must have at least one column")
  }

  cli::cli_alert_info("There are {ncol(sig_data)} signatures")

  sig_data <- as.data.frame(sig_data)
  sig_data[sig_data == "NA"] <- NA

  bb <- lapply(seq_len(ncol(sig_data)), function(i) {
    aa <- as.character(sig_data[, i])
    aa <- aa[!is.na(aa)]
    aa <- aa[aa != ""]
    unique(aa)
  })
  names(bb) <- names(sig_data)

  names(bb) <- gsub(names(bb), pattern = "\\ ", replacement = "_")

  my_signatures <- bb

  if (save_signature) {
    if (is.null(output_path)) {
      cli::cli_abort("{.arg output_path} must be provided when save_signature = TRUE")
    }
    save(my_signatures, file = paste0(output_path, ".RData"))
  }

  my_signatures
}
