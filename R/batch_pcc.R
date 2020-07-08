


#' Batch way to calculate the partial correlation coefficient
#'
#' @description  batch_pcc() provide a batch way to calculate the partial correlation coefficient between target gene and others when
#' controlling a third variable
#' @param data matrix;data signature matrix with multiple features
#' @param target a character vector; target name of group
#' @param feature character vectors; feature used to comparison
#' @param method vector; one of "pearson"(default), "spearman" or "kendall"
#' @param interference a numeric vector, matrix or data frame; variabes used to control when calculate the correlation between interesting variables.
#' @return
#' @export
#' @import purrr
#' @import ppcor
#' @import tidyverse
#' @author Rongfang Shen
#' @examples
#'
#' data("tcga_crc_exp")
#' data("tcga_crc_purity")
#' tcga_crc_exp <- t(tcga_crc_exp)
#' tcga_crc_exp <- data.frame(array = stringr：：str_sub(rownames(tcga_crc_exp), 1, 15), tcga_crc_exp)
#' dat <- merge(tcga_crc_purity, tcga_crc_exp, by = "array")
#' TP53cor <- batch_pcc(data = dat, target = "TP53",
#' feature = colnames(dat)[3:100], method = "pearson", interference = dat$purity)
batch_pcc <- function(data, target, feature, method = "pearson",
                      interference){
  feature <- feature[feature %in% colnames(data)];
  print(feature[1:2])
  aa <- data[, feature] %>% as.tibble() %>%
    map(pcor.test, y = data[,target], z = interference, method=method)
  pvalue <- aa %>% map_dbl("p.value")
  statistic <- aa %>% map_dbl("estimate")
  cc <- data.frame(sig_names = feature,
                 p.value = pvalue,
                 statistic = statistic)
  cc <- cc[order(cc$p.value, decreasing = F), ]
  cc$p.adj <- p.adjust(cc$p.value, method = "BH")
  cc$log10pvalue <- -1*log10(cc$p.value)
  rownames(cc) <- NULL
  cc$stars <- cut(cc$p.adj, breaks = c(-Inf,0.0001, 0.001, 0.01, 0.05,0.5, Inf),
                  label = c("****","***", "**", "*", "+",""))
  return(cc)

}
