#' Save Signature Data to File
#'
#' This function saves signature data to a specified file format, supporting CSV or RData.
#' It handles single signatures or lists of signatures, converting them to a data frame for storage.
#'
#' @param signatures Signature data: a list or single string of signatures.
#' @param format Output format: "csv" or "rdata". Default is "csv".
#' @param file.name Name of the output file without extension.
#'
#' @return A data frame containing the processed signature data, also saved to the specified file.
#' @export
#'
#' @examples
#' output_sig(signatures = signature_collection, format = "csv", file.name = "my_signatures")
output_sig <- function(signatures, format = "csv", file.name) {
  if (length(signatures) <= 1) {
    signatures <- as.data.frame(signatures)
    colnames(signatures) <- "signature genes"
  }

  if (length(signatures) >= 2) {
    signatures <- as.data.frame(t(do.call("rbind", signatures)))

    for (i in 1:ncol(signatures)) {
      index <- duplicated(signatures[, i])
      signatures[index, i] <- NA
    }
  }

  if (format == "csv") write.csv(signatures, file = paste0(file.name, ".csv"))

  if (format == "rdata") save(signatures, file = paste0(file.name, ".RData"))
  return(signatures)
}
