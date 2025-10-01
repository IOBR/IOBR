#' Batch Wilcoxon Rank-Sum Test Between Two Groups
#'
#' Performs Wilcoxon rank-sum tests on a dataset to compare the distribution of specified features between two groups. Computes p-values, adjusts for multiple testing, and ranks features by significance. Returns a data frame with feature names, p-values, adjusted p-values, log-transformed p-values, and significance stars.
#'
#' @param data Data frame. Input data for analysis.
#' @param target Character. Name of the column representing group labels. Default is "group".
#' @param feature Character vector or NULL. Names of features to analyze. If NULL, all continuous features are used. Default is NULL.
#' @param feature_manipulation Logical. Whether to apply custom feature manipulation. Default is FALSE.
#'
#' @return Data frame with statistical results for each feature (p-value, adjusted p-value, log p-value, significance stars).
#' @export
#' @import dplyr
#' @import tibble
#' @author Dongqiang Zeng
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Find microenvironmental scores associated with Gender
#' batch_wilcoxon(data = sig_stad, target = "Gender", feature = colnames(sig_stad)[69:ncol(sig_stad)])
batch_wilcoxon <- function(data, target = "group", feature = NULL, feature_manipulation = FALSE) {
  data <- as.data.frame(data)
  # change-name-of-group
  colnames(data)[which(colnames(data) == target)] <- "group"

  data <- data[!is.na(data$group), ]
  data <- data[!data$group == "", ]
  data$group <- as.character(data$group)
  group_names <- unique(data$group)

  if (is.null(feature)) {
    message(">>>-- `feature` must be specified, or all continuous features will be estimated...")

    index <- menu(c("all continuous features", "selected features "), title = " >>>-- Choose features:")
    if (index == 1) {
      feature <- colnames(data)
      feature <- feature[sapply(data, is.numeric)]
    } else {
      stop(">>>-- Please specify the features that you want to proceed...")
    }
  }
  feature <- feature[feature %in% colnames(data)]
  if (feature_manipulation) feature <- feature_manipulation(data = data, feature = feature, print_result = F)

  # if(!identical(group_names,c("High","Low"))) message(">>>--- `group_names` should be specified...")

  message(">>>-- Grouping information: ")
  print(table(data$group))
  ###########################################
  if (length(group_names) > 2) {
    print(table(data$group))
    stop("Variable has more than two levels...")
  }

  data <- data[, c("group", feature)]
  aa <- lapply(data[, feature], function(x) wilcox.test(x ~ data[, "group"], var.equal = F))
  result_mean <- data %>%
    dplyr::group_by(.$group) %>%
    dplyr::summarise_if(is.numeric, mean)

  rownames(result_mean) <- NULL
  result_mean <- result_mean %>%
    tibble::column_to_rownames(., var = ".$group") %>%
    base::as.data.frame() %>%
    t() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(., var = "sig_names")

  result_mean <- as.data.frame(result_mean)
  group_names <- group_names[order(group_names)]
  colnames(result_mean)[2:3] <- group_names
  result_mean$statistic <- result_mean[, 2] - result_mean[, 3]

  cc <- data.frame(
    sig_names = feature,
    p.value = sapply(aa, getElement, name = "p.value")
  )

  cc <- cc %>%
    full_join(result_mean, by = "sig_names") %>%
    dplyr::arrange(p.value) %>%
    dplyr::mutate(p.adj = p.adjust(.$p.value, method = "BH")) %>%
    dplyr::mutate(log10pvalue = log10(.$p.value) * -1) %>%
    dplyr::mutate(stars = cut(.$p.value,
      breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.5, Inf),
      label = c("****", "***", "**", "*", "+", "")
    ))
  cc <- tibble::as_tibble(cc)
  return(cc)
}
