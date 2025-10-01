#' Identify High-Variance Features from Statistical Results
#'
#' Selects top variable (up- and down-regulated) features based on adjusted p-value and log fold-change thresholds.
#'
#' @param result Data frame or tibble. Statistical results containing feature, adjusted p-value, and logFC columns.
#' @param target Character. Column name of feature identifiers.
#' @param name_padj Character. Adjusted p-value column name. Default "padj".
#' @param name_logfc Character. log2 fold-change column name.
#' @param n Integer. Number of top up and top down features to select. Default 10.
#' @param padj_cutoff Numeric. Adjusted p-value threshold. Default 1.
#' @param data_type Character or NULL. If "survival", adjusts logFC interpretation.
#' @param logfc_cutoff Numeric. Absolute log2 fold-change threshold. Default 0.
#'
#' @return Character vector of selected feature names (combined up and down sets).
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' result_data <- data.frame(gene = c("Gene1","Gene2","Gene3"), padj = c(0.01,0.2,0.05), logfc = c(-2,1.5,-3))
#' high_var_fea(result = result_data, target = "gene", name_padj = "padj", name_logfc = "logfc", n = 2,
#'              padj_cutoff = 0.05, logfc_cutoff = 1.5)
high_var_fea <- function(result, target, name_padj = "padj", padj_cutoff = 1, name_logfc, logfc_cutoff = 0, n = 10, data_type = NULL) {
  colnames(result)[which(colnames(result) == name_padj)] <- "padj"
  colnames(result)[which(colnames(result) == name_logfc)] <- "logfc"
  #############################################
  if (!is.null(data_type)) {
    if (data_type == "survival") result[, "logfc"] <- result[, "logfc"] - 1
  }
  #############################################
  colnames(result)[which(colnames(result) == target)] <- "target"

  topVarFeature1 <- result %>%
    # filter(!is.na(target)) %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::filter(logfc < -abs(logfc_cutoff)) %>%
    dplyr::select(target, padj, logfc)

  topVarFeature2 <- result %>%
    # filter(!is.na(target)) %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::filter(logfc > abs(logfc_cutoff)) %>%
    dplyr::select(target, padj, logfc)

  if (dim(topVarFeature1)[1] < n) {
    message(paste0(">>> The cutoff was too strict for down regulated features, only ", dim(topVarFeature1)[1], "  variables were found in the result"))
  }

  if (dim(topVarFeature2)[1] < n) {
    message(paste0(">>> The cutoff was too strict for up regulated features, only ", dim(topVarFeature2)[1], "  variables were found in the result"))
  }

  if (dim(topVarFeature1)[1] == 0) message(paste0(">>> No down regulated features was found"))

  if (dim(topVarFeature2)[1] == 0) message(paste0(">>> No up regulated features was found"))

  topVarFeature <- c(
    as.character(topVarFeature1$target)[1:n],
    as.character(topVarFeature2$target)[1:n]
  )
  topVarFeature <- topVarFeature[!is.na(topVarFeature)]
  return(topVarFeature)
}
