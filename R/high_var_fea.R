#' Identify High-Variance Features from Statistical Results
#'
#' @description
#' Selects top variable (up- and down-regulated) features based on adjusted
#' p-value and log fold-change thresholds.
#'
#' @param result Data frame or tibble. Statistical results containing feature,
#'   adjusted p-value, and logFC columns.
#' @param target Character. Column name of feature identifiers.
#' @param name_padj Character. Adjusted p-value column name. Default is `"padj"`.
#' @param padj_cutoff Numeric. Adjusted p-value threshold. Default is 1.
#' @param name_logfc Character. log2 fold-change column name.
#' @param logfc_cutoff Numeric. Absolute log2 fold-change threshold. Default is 0.
#' @param n Integer. Number of top up and top down features to select.
#'   Default is 10.
#' @param data_type Character or `NULL`. If `"survival"`, adjusts logFC
#'   interpretation. Default is `NULL`.
#'
#' @return Character vector of selected feature names (combined up and down sets).
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' result_data <- data.frame(
#'   gene = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
#'   padj = c(0.01, 0.02, 0.05, 0.001, 0.03),
#'   logfc = c(-2, 1.5, -3, 2.5, 0.5)
#' )
#' high_var_fea(
#'   result = result_data,
#'   target = "gene",
#'   name_padj = "padj",
#'   name_logfc = "logfc",
#'   n = 2,
#'   padj_cutoff = 0.05,
#'   logfc_cutoff = 1.5
#' )
high_var_fea <- function(result, target, name_padj = "padj", padj_cutoff = 1,
                         name_logfc, logfc_cutoff = 0, n = 10, data_type = NULL) {
  if (!is.data.frame(result)) {
    cli::cli_abort("{.arg result} must be a data frame")
  }
  if (!target %in% colnames(result)) {
    cli::cli_abort("Column {.val {target}} not found in {.arg result}")
  }
  if (!name_padj %in% colnames(result)) {
    cli::cli_abort("Column {.val {name_padj}} not found in {.arg result}")
  }
  if (!name_logfc %in% colnames(result)) {
    cli::cli_abort("Column {.val {name_logfc}} not found in {.arg result}")
  }

  colnames(result)[colnames(result) == name_padj] <- "padj"
  colnames(result)[colnames(result) == name_logfc] <- "logfc"
  colnames(result)[colnames(result) == target] <- "target"

  if (!is.null(data_type) && data_type == "survival") {
    result[, "logfc"] <- result[, "logfc"] - 1
  }

  topVarFeature1 <- result %>%
    dplyr::filter(.data$padj < padj_cutoff) %>%
    dplyr::arrange(.data$padj) %>%
    dplyr::filter(.data$logfc < -abs(logfc_cutoff)) %>%
    dplyr::select(.data$target, .data$padj, .data$logfc)

  topVarFeature2 <- result %>%
    dplyr::filter(.data$padj < padj_cutoff) %>%
    dplyr::arrange(.data$padj) %>%
    dplyr::filter(.data$logfc > abs(logfc_cutoff)) %>%
    dplyr::select(.data$target, .data$padj, .data$logfc)

  if (nrow(topVarFeature1) < n) {
    cli::cli_alert_warning(
      "Cutoff too strict for down-regulated features, only {nrow(topVarFeature1)} found"
    )
  }
  if (nrow(topVarFeature2) < n) {
    cli::cli_alert_warning(
      "Cutoff too strict for up-regulated features, only {nrow(topVarFeature2)} found"
    )
  }
  if (nrow(topVarFeature1) == 0) {
    cli::cli_alert_warning("No down-regulated features found")
  }
  if (nrow(topVarFeature2) == 0) {
    cli::cli_alert_warning("No up-regulated features found")
  }

  topVarFeature <- c(
    as.character(topVarFeature1$target)[seq_len(min(n, nrow(topVarFeature1)))],
    as.character(topVarFeature2$target)[seq_len(min(n, nrow(topVarFeature2)))]
  )
  topVarFeature <- topVarFeature[!is.na(topVarFeature)]

  topVarFeature
}
