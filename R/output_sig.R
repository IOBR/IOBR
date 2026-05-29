#' Save Signature Data to File
#'
#' @description
#' Saves signature data to a specified file format, supporting CSV or RData.
#' Handles single signatures or lists of signatures, converting them to a
#' data frame for storage.
#'
#' @param signatures Signature data: a list or single string of signatures.
#' @param format Character. Output format: "csv" or "rdata". Default is "csv".
#' @param file.name Character. Name of the output file without extension.
#'
#' @return Data frame containing the processed signature data, also saved to
#'   the specified file.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' # Simulate signatures
#' set.seed(123)
#' sim_sigs <- list(
#'   Signature1 = paste0("Gene", 1:50),
#'   Signature2 = paste0("Gene", 51:100),
#'   Signature3 = paste0("Gene", 101:150)
#' )
#' tmpfile <- tempfile(fileext = ".csv")
#' output_sig(
#'   signatures = sim_sigs, format = "csv",
#'   file.name = tools::file_path_sans_ext(tmpfile)
#' )
output_sig <- function(signatures, format = c("csv", "rdata"), file.name) {
  format <- rlang::arg_match(format)

  if (missing(file.name)) {
    cli::cli_abort("{.arg file.name} must be provided.")
  }

  if (is.null(signatures) || length(signatures) == 0) {
    cli::cli_abort("{.arg signatures} cannot be NULL or empty.")
  }

  if (length(signatures) == 1) {
    signatures <- as.data.frame(signatures)
    colnames(signatures) <- "signature genes"
  } else {
    signatures <- as.data.frame(t(do.call("rbind", signatures)))

    for (i in seq_len(ncol(signatures))) {
      signatures[duplicated(signatures[, i]), i] <- NA
    }
  }

  if (format == "csv") {
    write.csv(signatures, file = paste0(file.name, ".csv"), row.names = FALSE)
    cli::cli_alert_success("Signature data saved to {.file {file.name}.csv}")
  } else {
    save(signatures, file = paste0(file.name, ".RData"))
    cli::cli_alert_success("Signature data saved to {.file {file.name}.RData}")
  }

  signatures
}
