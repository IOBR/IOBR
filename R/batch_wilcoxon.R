


#' Batch to conduct wilcoxon test with two groups
#'
#' @param data signature matrix with two groups
#' @param target name of group: with `high` and `low` category
#' @param feature feature used to comparison
#' @return
#' @export
#' @import dplyr
#' @import tibble
#' @author Dongqiang Zeng
#' @examples
batch_wilcoxon<-function(data = data,target = "group", feature = feature){

  feature<-feature[feature%in%colnames(data)];print(feature[1:2])

  #change-name-of-group
  colnames(data)[which(colnames(data)==target)]<-"group"
  aa<-lapply(data[,feature], function(x) wilcox.test(x ~ data[,"group"],var.equal = F))

  result_mean<-data %>% group_by(.$group) %>%
    summarise_if(is.numeric,mean) %>%
    column_to_rownames(.,var = ".$group") %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(.,var = "sig_names") %>%
    mutate(statistic = .$High - .$Low)

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

