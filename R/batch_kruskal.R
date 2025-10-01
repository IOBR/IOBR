#' Batch Kruskal-Wallis Test
#'
#' @description
#' Performs Kruskal-Wallis rank sum tests on multiple continuous features across
#' different groups. Computes p-values, adjusts for multiple testing, and ranks
#' features by significance.
#'
#' @param data Data frame containing the dataset for analysis.
#' @param group Character string specifying the name of the grouping variable.
#' @param feature Character vector specifying the names of feature variables to test.
#'   If \code{NULL}, the user is prompted to select all continuous features or specify
#'   features manually. Default is \code{NULL}.
#' @param feature_manipulation Logical indicating whether to apply feature manipulation
#'   to filter valid features. Default is \code{FALSE}.
#'
#' @return Tibble containing the following columns for each feature:
#' \itemize{
#'   \item \code{sig_names}: Feature name
#'   \item \code{p.value}: Raw p-value from Kruskal-Wallis test
#'   \item \code{p.adj}: Adjusted p-value (Benjamini-Hochberg method)
#'   \item \code{log10pvalue}: Negative log10-transformed p-value
#'   \item \code{stars}: Significance stars based on p-value thresholds
#' }
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Test features associated with TCGA molecular subtype
#' batch_kruskal(data = sig_stad, group = "Subtype",
#'               feature = colnames(sig_stad)[69:ncol(sig_stad)])
batch_kruskal <- function(data, group, feature = NULL, feature_manipulation = FALSE) {
  data <- as.data.frame(data)

  feature <- feature[feature %in% colnames(data)]
  data <- data[, c(group, feature)]

  # change-name-of-group
  colnames(data)[which(colnames(data) == group)] <- "group"

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

  if (feature_manipulation) feature <- feature_manipulation(data = data, feature = feature, print_result = F)

  data <- data[, c("group", feature)]
  # if(!identical(group_names,c("High","Low"))) message(">>>--- `group_names` should be specified...")

  message(">>>-- Grouping information: ")
  print(table(data$group))
  ###########################################

  if (length(group_names) < 3) {
    print(table(data$group))
    stop("Variable has less than three levels...")
  }


  aa <- lapply(data[, feature], function(x) kruskal.test(x ~ data[, "group"]))
  res <- data.frame(
    p.value = sapply(aa, getElement, name = "p.value"),
    sig_names = feature,
    statistic = sapply(aa, getElement, name = "statistic")
  )
  res$p.adj <- p.adjust(res$p.value, method = "BH", n = length(res$p.value))
  res <- res[order(res$p.adj, decreasing = F), ]

  # writexl::write_xlsx(res, paste0(file_name$abspath, prefix, "0-statistical-res-with-",group,".xlsx"))
  #######################################

  result_mean <- data %>%
    dplyr::group_by(.$group) %>%
    dplyr::summarise_if(is.numeric, mean)

  # mean of each group
  rownames(result_mean) <- NULL
  result_mean <- result_mean %>%
    tibble::column_to_rownames(., var = ".$group") %>%
    base::as.data.frame() %>%
    t() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(., var = "sig_names")

  result_mean <- as.data.frame(result_mean)
  group_names <- group_names[order(group_names)]
  colnames(result_mean)[2:ncol(result_mean)] <- group_names

  result_mean$mean <- rowSums(result_mean[, 2:ncol(result_mean)]) / length(group_names)
  result_mean[, group_names] <- apply(result_mean[, group_names], 2, function(x) x - result_mean$mean)
  # result_mean$statistic<- result_mean[,2] - result_mean[,3]

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

  # input<-extract_sc_data(sce = sces_merged, vars = feas, slot = "data", assay = "RNA")
  # source("E:/18-Github/Organization/IOBR/R/batch_kruskal.R")
  # batch_kruskal(data = input, group = "celltype_modified", feature = feas, feature_manipulation = FALSE)
}
