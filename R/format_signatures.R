#' Transform Signature Data into List Format
#'
#' This function converts signature data from a data frame (with signatures as columns and genes as rows)
#' into a list format suitable for IOBR functions, replacing NA values appropriately.
#'
#' @param sig_data Data frame with signature names as columns and genes in rows; use NA for missing values.
#' @param save_signature Logical indicating whether to save the signature list as RData. Default is FALSE.
#' @param output_name String for the output RData file name. Default is "signatures".
#'
#' @return A list of signatures.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Load signature data frame
#' data("sig_excel", package = "IOBR")
#' # Transform into gene list for IOBR functions
#' format_signatures(sig_excel)
format_signatures <- function(sig_data, save_signature = FALSE, output_name = "signatures") {
  message(paste0(">>> There are ", dim(sig_data)[2], "  signatures >>>"))

  sig_data <- as.data.frame(sig_data)
  sig_data[sig_data == "NA"] <- NA
  #' reading each column and transfer it to a list
  bb <- as.list(NULL)
  for (i in 1:ncol(sig_data)) {
    aa <- as.character(sig_data[, i])
    aa <- list(aa)
    names(aa) <- names(sig_data[i])
    bb <- append(bb, aa)
  }
  bb <- lapply(bb, function(x) na.omit(x))
  bb <- lapply(bb, function(x) as.character(x))
  bb <- lapply(bb, function(x) unique(x))
  bb <- lapply(bb, function(x) x[!x == ""])
  #' standerdized the name of list
  # names(bb)<-gsub(names(bb),pattern = "\\.",replacement = "_")
  names(bb) <- gsub(names(bb), pattern = "\\ ", replacement = "_")
  # names(bb)<-gsub(names(bb),pattern = "\\-",replacement = "_")
  # ########################################
  my_signatures <- bb
  ##########################################
  if (save_signature == TRUE) {
    save(my_signatures, file = paste0(output_name, ".RData"))
  }
  return(my_signatures)
}
