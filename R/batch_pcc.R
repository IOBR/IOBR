


#' Batch way to calculate the partial correlation coefficient
#'
#' @description  batch_pcc() provide a batch way to calculate the partial correlation coefficient between target gene and others when
#' controlling a third variable
#' @param pdata_group matrix;data signature matrix with multiple features
#' @param feature_data
#' @param id1
#' @param id2
#' @param interferenceid character vectors; vector used to control
#' @param target character vectors; target name of group
#' @param method vector; one of "pearson"(default), "spearman" or "kendall"
#' @return
#' @export
#' @author Rongfang Shen
#' @examples
#' pdata_group <- imvigor210_pdata[, c("ID", "TumorPurity", "Pan_F_TBRs")] %>%
#' rename(target = Pan_F_TBRs) %>% mutate(target = as.numeric(target))
#' res <- batch_pcc(pdata_group = pdata_group, id1 = "ID", feature_data = imvigor210_sig,
#' id2 = "ID", interferenceid = "TumorPurity",
#' target = "target", method = "pearson")
batch_pcc <- function(pdata_group, id1 = "ID",
                      feature_data, id2 = "ID",
                      interferenceid,
                      target, method = "pearson"){
  dat <- merge(pdata_group, feature_data, by.x = id1, by.y = id2)
  features <- setdiff(colnames(feature_data), id2)
  aa <- dat[, features] %>%as_tibble() %>%
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

