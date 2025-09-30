#' Extract most variable features form statistical results
#' @description This function is designed to identify high variability features (genes or markers) based on certain criteria from a given dataset. It takes as input the dataset and several parameters, such as the significance cutoff for adjusted p-values, the fold change cutoff for log-fold changes, the number of top variables to select, and the data type. The function then filters the dataset based on these criteria and selects the top variable features that meet the specified conditions.
#' @param result a tibble or data frame
#' @param target The column name of the target variable in the dataset.
#' @param name_padj The column name representing the adjusted p-value in the dataset. Default is "padj".
#' @param name_logfc he column name representing the log-fold change in the dataset. Default is NULL.
#' @param n The number of top variable features to select. Default is 10.
#' @param padj_cutoff The significance cutoff for adjusted p-values. Values below this threshold are considered significant. Default is 1.
#' @param data_type The type of data being analyzed. Options include "survival" and other relevant data types. Default is NULL.
#' @param logfc_cutoff The fold change cutoff for log-fold changes. Values above or below this threshold are considered significant. Default is 0.
#'
#' @author Dongqiang Zeng
#' @return A vector of feature names identified as the most variable, combining both upregulated and downregulated features.
#' @export
#'
#' @examples
#' # Assume 'result_data' is a data frame with statistical analysis results
#' result_data <- data.frame(
#'   gene = c("Gene1", "Gene2", "Gene3"),
#'   padj = c(0.01, 0.2, 0.05),
#'   logfc = c(-2, 1.5, -3)
#' )
#' high_var_features <- high_var_fea(
#'   result = result_data,
#'   target = "gene",
#'   name_padj = "padj",
#'   name_logfc = "logfc",
#'   n = 2,
#'   padj_cutoff = 0.05,
#'   logfc_cutoff = 1.5
#' )
#' print(high_var_features)
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
