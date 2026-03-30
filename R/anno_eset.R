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
#' \enumerate{
#'   \item Filters probes with missing symbols or labeled as `"NA_NA"`
#'   \item Matches probes between expression set and annotation data
#'   \item Merges annotation with expression data
#'   \item Handles duplicate gene symbols using specified aggregation method
#'   \item Removes rows with all zeros, all NAs, or missing values in the first column
#' }
#'
#' @param eset Expression matrix or ExpressionSet object containing gene
#'   expression data.
#' @param annotation Data frame containing annotation information for probes.
#'   Built-in options include `anno_hug133plus2`, `anno_rnaseq`, and
#'   `anno_illumina`.
#' @param symbol Character string specifying the column name in `annotation`
#'   that represents gene symbols. Default is `"symbol"`.
#' @param probe Character string specifying the column name in `annotation`
#'   that represents probe identifiers. Default is `"probe_id"`.
#' @param method Character string specifying the aggregation method for
#'   duplicate gene symbols. Options are `"mean"`, `"sum"`, or `"sd"`.
#'   Default is `"mean"`.
#'
#' @return Annotated and cleaned gene expression matrix with gene symbols as
#'   row names.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' # Annotate Affymetrix microarray data
#' eset_gse62254 <- load_data("eset_gse62254")
#' anno_hug133plus2 <- load_data("anno_hug133plus2")
#' eset <- anno_eset(eset = eset_gse62254, annotation = anno_hug133plus2)
#' head(eset)
#' }
#'
#' \donttest{
#' # Annotate RNA-seq data with Ensembl IDs
#' eset_stad <- load_data("eset_stad")
#' anno_grch38 <- load_data("anno_grch38")
#' eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")
#' head(eset)
#' }
anno_eset <- function(eset,
                      annotation,
                      symbol = "symbol",
                      probe = "probe_id",
                      method = "mean") {

  # Input validation
  if (is.null(eset)) {
    cli::cli_abort("{.arg eset} cannot be NULL.")
  }
  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort(c(
      "Invalid type for {.arg eset}.",
      "i" = "Expected a matrix or data frame, got {.cls {class(eset)}}."
    ))
  }
  if (nrow(eset) == 0) {
    cli::cli_abort("{.arg eset} has no rows.")
  }
  if (is.null(annotation) || !is.data.frame(annotation)) {
    cli::cli_abort(c(
      "Invalid type for {.arg annotation}.",
      "i" = "Expected a data frame, got {.cls {class(annotation)}}."
    ))
  }
  if (!symbol %in% colnames(annotation)) {
    cli::cli_abort("Symbol column {.val {symbol}} not found in annotation.")
  }
  if (!probe %in% colnames(annotation)) {
    cli::cli_abort("Probe column {.val {probe}} not found in annotation.")
  }
  method <- rlang::arg_match(method, c("mean", "sum", "sd"))

  # Prepare annotation
  annotation <- as.data.frame(annotation)
  colnames(annotation)[colnames(annotation) == symbol] <- "symbol"
  colnames(annotation)[colnames(annotation) == probe] <- "probe_id"
  annotation <- annotation[, c("probe_id", "symbol")]

  # Filter out NA and invalid symbols
  valid_symbols <- !is.na(annotation$symbol) & annotation$symbol != "NA_NA"
  annotation <- annotation[valid_symbols, , drop = FALSE]

  cli::cli_alert_info("Row number of original eset: {nrow(eset)}")

  # Match probes
  matched_probes <- rownames(eset) %in% annotation$probe_id
  anno_count <- sum(matched_probes) / length(matched_probes)
  cli::cli_alert_success(
    "{format(100 * anno_count, digits = 2)}% of probes in expression set were annotated"
  )

  if (sum(matched_probes) == 0) {
    cli::cli_abort(c(
      "No probes matched between eset and annotation.",
      "i" = "Check that probe ID formats match between eset rownames and annotation.",
      "*" = "Example eset rowname: {.val {rownames(eset)[1]}}",
      "*" = "Example annotation probe: {.val {annotation$probe_id[1]}}"
    ))
  }

  # Filter eset to matched probes
  eset <- eset[matched_probes, , drop = FALSE]
  annotation <- annotation[annotation$probe_id %in% rownames(eset), , drop = FALSE]

  # Merge
  eset <- as.data.frame(eset)
  eset <- tibble::rownames_to_column(eset, var = "id")
  eset <- merge(annotation, eset, by.x = "probe_id", by.y = "id", all = FALSE)
  eset$probe_id <- NULL

  # Handle duplicates
  n_dups <- nrow(eset) - length(unique(eset$symbol))
  if (n_dups > 0) {
    cli::cli_alert_info(
      "Found {n_dups} duplicate symbol{?s}, using {.val {method}} method"
    )

    agg_col <- setdiff(colnames(eset), "symbol")

    order_index <- switch(method,
      "mean" = rowMeans(eset[, agg_col, drop = FALSE], na.rm = TRUE),
      "sum"  = rowSums(eset[, agg_col, drop = FALSE], na.rm = TRUE),
      "sd"   = apply(eset[, agg_col, drop = FALSE], 1, stats::sd, na.rm = TRUE)
    )

    eset <- eset[order(order_index, decreasing = TRUE), , drop = FALSE]
    eset <- dplyr::distinct(eset, symbol, .keep_all = TRUE)
  }

  eset <- tibble::column_to_rownames(eset, var = "symbol")

  # Remove uninformative rows
  eset <- eset[rowSums(eset == 0) < ncol(eset), , drop = FALSE]
  eset <- eset[rowSums(is.na(eset)) < ncol(eset), , drop = FALSE]
  eset <- eset[!is.na(eset[, 1]), , drop = FALSE]

  cli::cli_alert_info(
    "Row number after filtering duplicated gene symbol: {nrow(eset)}"
  )

  as.matrix(eset)
}
