#' Merge Expression Sets by Row Names
#'
#' @description
#' Merges two or three expression sets (matrices or data frames) by row names
#' (gene symbols), removing duplicates. The function ensures common genes
#' across all input expression sets are retained.
#'
#' @param eset1 First expression set (matrix or data frame with row names).
#' @param eset2 Second expression set (matrix or data frame with row names).
#' @param eset3 Optional third expression set. Default is `NULL`.
#'
#' @return Merged expression set (data frame) with duplicates removed.
#'   Row names correspond to gene symbols.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' # Load example data
#' eset_stad <- load_data("eset_stad")
#'
#' # Create mock expression sets with common genes
#' common_genes <- c("TP53", "BRCA1", "EGFR", "MYC")
#' eset1 <- matrix(rnorm(12),
#'   nrow = 4,
#'   dimnames = list(common_genes, paste0("S", 1:3))
#' )
#' eset2 <- matrix(rnorm(16),
#'   nrow = 4,
#'   dimnames = list(common_genes, paste0("S", 4:7))
#' )
#'
#' # Merge two expression sets
#' merged_eset <- merge_eset(eset1, eset2)
#' print(dim(merged_eset))
#' }
merge_eset <- function(eset1, eset2, eset3 = NULL) {
  # Input validation
  if (!is.matrix(eset1) && !is.data.frame(eset1)) {
    cli::cli_abort(c(
      "Invalid type for {.arg eset1}.",
      "i" = "Expected a matrix or data frame, got {.cls {class(eset1)}}."
    ))
  }
  if (!is.matrix(eset2) && !is.data.frame(eset2)) {
    cli::cli_abort(c(
      "Invalid type for {.arg eset2}.",
      "i" = "Expected a matrix or data frame, got {.cls {class(eset2)}}."
    ))
  }

  # Convert to data frames and ensure row names
  eset1 <- as.data.frame(eset1)
  eset2 <- as.data.frame(eset2)

  if (is.null(rownames(eset1)) || is.null(rownames(eset2))) {
    cli::cli_abort("Input expression sets must have row names (gene symbols).")
  }

  # First merge
  common_genes_12 <- intersect(rownames(eset1), rownames(eset2))
  cli::cli_alert_info(
    "Common genes between eset1 and eset2: {length(common_genes_12)}"
  )

  if (length(common_genes_12) == 0) {
    cli::cli_abort("No common genes found between eset1 and eset2.")
  }

  eset <- merge(eset1, eset2, by = "row.names", all = FALSE)
  eset <- remove_duplicate_genes(eset = eset, column_of_symbol = "Row.names")

  # Optional third merge
  if (!is.null(eset3)) {
    if (!is.matrix(eset3) && !is.data.frame(eset3)) {
      cli::cli_abort(c(
        "Invalid type for {.arg eset3}.",
        "i" = "Expected a matrix or data frame, got {.cls {class(eset3)}}."
      ))
    }
    eset3 <- as.data.frame(eset3)

    if (is.null(rownames(eset3))) {
      cli::cli_abort("eset3 must have row names (gene symbols).")
    }

    common_genes_3 <- intersect(rownames(eset), rownames(eset3))
    cli::cli_alert_info(
      "Common genes after merging with eset3: {length(common_genes_3)}"
    )

    if (length(common_genes_3) == 0) {
      cli::cli_abort("No common genes found with eset3.")
    }

    eset <- merge(eset, eset3, by = "row.names", all = FALSE)
    eset <- remove_duplicate_genes(eset = eset, column_of_symbol = "Row.names")
  }

  cli::cli_alert_success(
    "Final merged expression set: {nrow(eset)} genes x {ncol(eset)} samples"
  )

  eset
}
