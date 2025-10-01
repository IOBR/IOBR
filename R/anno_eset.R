#' Annotate Gene Expression Matrix and Remove Duplicated Genes
#'
#' @description
#' Annotates an expression matrix with gene symbols using provided annotation data,
#' filters out missing or invalid symbols, handles duplicate gene entries, and
#' removes uninformative rows. The function supports multiple aggregation methods
#' for resolving duplicate gene symbols.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Filters probes with missing symbols or labeled as "NA_NA"
#'   \item Matches probes between expression set and annotation data
#'   \item Merges annotation with expression data
#'   \item Handles duplicate gene symbols using specified aggregation method
#'   \item Removes rows with all zeros, all NAs, or missing values in the first column
#' }
#'
#' @param eset Expression matrix or ExpressionSet object containing gene expression data.
#' @param annotation Data frame containing annotation information for probes. Built-in
#'   options include \code{anno_hug133plus2}, \code{anno_rnaseq}, and \code{anno_illumina}.
#' @param symbol Character string specifying the column name in \code{annotation} that
#'   represents gene symbols. Default is \code{"symbol"}.
#' @param probe Character string specifying the column name in \code{annotation} that
#'   represents probe identifiers. Default is \code{"probe_id"}.
#' @param method Character string specifying the aggregation method for duplicate gene
#'   symbols. Options are \code{"mean"}, \code{"sum"}, or \code{"sd"}. Default is
#'   \code{"mean"}.
#'
#' @return Annotated and cleaned gene expression matrix with gene symbols as row names.
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Annotate Affymetrix microarray data
#' data(eset_gse62254, package = "IOBR")
#' eset <- anno_eset(eset = eset_gse62254, annotation = anno_hug133plus2)
#'
#' # Annotate RNA-seq data with Ensembl IDs
#' data(eset_stad, package = "IOBR")
#' eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
anno_eset <- function(eset, annotation, symbol = "symbol", probe = "probe_id", method = "mean") {
  annotation <- as.data.frame(annotation)
  colnames(annotation)[which(colnames(annotation) == symbol)] <- "symbol"
  colnames(annotation)[which(colnames(annotation) == probe)] <- "probe_id"
  annotation <- annotation[, c("probe_id", "symbol")]
  annotation <- annotation[!annotation$symbol == "NA_NA", ]
  annotation <- annotation[!is.na(annotation$symbol), ]

  message(paste0("Row number of original eset: "))
  message(paste0(">>>>  ", dim(eset)[1]))

  annotation <- annotation[!is.na(annotation$symbol), ]

  anno_count <- length(rownames(eset)[rownames(eset) %in% annotation$probe_id]) / length(rownames(eset))
  message(paste0(paste0(sprintf(">>> %1.2f%%", 100 * anno_count), " of probe in expression set was annotated")))

  annotation <- annotation[annotation$probe_id %in% rownames(eset), ]
  eset <- eset[rownames(eset) %in% annotation$probe_id, ]

  eset <- as.data.frame(eset)
  eset <- tibble::rownames_to_column(eset, var = "id")
  eset <- merge(annotation, eset, by.x = "probe_id", by.y = "id", all = F)
  eset <- eset[, -1]

  eset <- as.data.frame(eset)
  rownames(eset) <- NULL

  symbol <- "symbol"
  dups <- dim(eset)[1] - length(unique(eset[, symbol]))

  if (dups == 0) {
    eset <- tibble::column_to_rownames(eset, var = symbol)
  } else {
    if (method == "mean") {
      order_index <- apply(eset[, setdiff(colnames(eset), symbol)], 1, function(x) mean(x, na.rm = T))
      eset <- eset[order(order_index, decreasing = T), ]
      eset <- eset %>%
        dplyr::distinct(!!sym(symbol), .keep_all = TRUE) %>%
        tibble::column_to_rownames(., var = symbol)
    } else if (method == "sd") {
      order_index <- apply(eset[, setdiff(colnames(eset), symbol)], 1, function(x) sd(x, na.rm = T))
      eset <- eset[order(order_index, decreasing = T), ]
      eset <- eset %>%
        distinct(!!sym(symbol), .keep_all = TRUE) %>%
        tibble::column_to_rownames(., var = symbol)
    } else if (method == "sum") {
      order_index <- apply(eset[, setdiff(colnames(eset), symbol)], 1, function(x) sum(x, na.rm = T))
      eset <- eset[order(order_index, decreasing = T), ]
      eset <- eset %>%
        distinct(!!sym(symbol), .keep_all = TRUE) %>%
        tibble::column_to_rownames(., var = symbol)
    }
  }

  quit <- rowSums(eset == 0) == ncol(eset)
  eset <- eset[!quit, ]

  quit <- rowSums(is.na(eset)) == ncol(eset)
  eset <- eset[!quit, ]

  quit <- is.na(eset[, 1])
  eset <- eset[!quit, ]

  message(paste0("Row number after filtering duplicated gene symbol: "))
  message(paste0(">>>>  ", dim(eset)[1]))

  return(eset)
}
