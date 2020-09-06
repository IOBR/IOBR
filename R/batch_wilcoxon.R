


#' Batch to conduct wilcoxon test with two groups
#'
#' @param data signature matrix with two groups
#' @param target name of group:
#' @param group_names  with `high` and `low` default
#' @param feature feature used to comparison
#'
#' @return
#' @export
#' @import dplyr
#' @import tibble
#' @author Dongqiang Zeng
#' @examples
batch_wilcoxon<-function(data = data,target = "group",group_names = c("High","Low"), feature = feature){

  feature<-feature[feature%in%colnames(data)]

  #change-name-of-group
  colnames(data)[which(colnames(data)==target)]<-"group"
  aa<-lapply(data[,feature], function(x) wilcox.test(x ~ data[,"group"],var.equal = F))

  result_mean<-data %>% group_by(.$group) %>%
    summarise_if(is.numeric,mean) %>%
    column_to_rownames(.,var = ".$group") %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(.,var = "sig_names")

  result_mean<-as.data.frame(result_mean)
  group_names<-group_names[order(group_names)]
  colnames(result_mean)[2:3]<-group_names
  result_mean$statistic<- result_mean[,2] - result_mean[,3]

  cc<-data.frame(sig_names=feature,
                 p.value = sapply(aa, getElement, name = "p.value"))

  cc<-cc %>%  full_join(result_mean,by = "sig_names") %>%
    arrange(p.value) %>%
    mutate(p.adj = p.adjust(.$p.value,method = "BH") ) %>%
    mutate(log10pvalue = log10(.$p.value)* -1) %>%
    mutate(stars = cut(.$p.value, breaks=c(-Inf,0.0001, 0.001, 0.01, 0.05,0.5, Inf),
                       label=c("****","***", "**", "*", "+","")))
  return(cc)
}

