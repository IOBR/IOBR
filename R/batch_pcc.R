

#' Batch way to calculate the partial correlation coefficient
#'
#' @description  batch_pcc() provide a batch way to calculate the partial correlation coefficient between feature and others when
#' controlling a third variable
#' @param pdata_group matrix;data signature matrix with multiple features

#' @param input A data frame containing both feature variables and the interference variable.
#' @param interferenceid The name of the column in the feature_data data frame representing the interference variable.
#' @param target The name of the column in the input data frame representing the target variable for correlation.
#' @param features A character vector specifying the names of the feature variables.
#' @param method The correlation method to be used. Default value is "pearson"; options are "pearson", "spearman", or "kendall".
#' 
#' @return A tibble containing the feature names, partial correlation coefficients, p-values, adjusted p-values, log10 p-values, and significance stars.
#' @export
#' @author Rongfang Shen
#' @examples
#' # Loading TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Finding Pan_F_TBRs associated signature score excluding the effects of tumour purity.
#' res <- batch_pcc(input = sig_stad, interferenceid = "TumorPurity_estimate", target = "Pan_F_TBRs", method = "pearson", features = colnames(sig_stad)[70:ncol(sig_stad)])
batch_pcc <- function(input, interferenceid, target, features, method = "pearson"){

  dat <- input
  features <- setdiff(features, c(unique(interferenceid, target)))

  aa <- dat[, features] %>% tibble::as_tibble() %>%
    map(ppcor::pcor.test, y = dat[,target], z = dat[, interferenceid], method=method)
  pvalue <- aa %>% purrr::map_dbl("p.value")
  statistic <- aa %>% purrr::map_dbl("estimate")
  cc <- data.frame(sig_names = features,
                   p.value = pvalue,
                   statistic = statistic)
  cc <- cc[order(cc$p.value, decreasing = F), ]
  cc$p.adj <- p.adjust(cc$p.value, method = "BH")
  cc$log10pvalue <- -1*log10(cc$p.value)
  rownames(cc) <- NULL
  cc$stars <- cut(cc$p.adj, breaks = c(-Inf,0.0001, 0.001, 0.01, 0.05,0.5, Inf),
                  label = c("****","***", "**", "*", "+",""))
  cc<-tibble::as_tibble(cc)
  return(cc)
}
