
#' Batch to conduct wilcoxon test between two groups
#'
#' @description The batch_wilcoxon function performs Wilcoxon rank-sum tests on a given data set to compare the distribution of a specified feature between two groups. It computes the p-values and ranks the significant features based on the p-values. The function returns a data frame with the feature names, p-values, adjusted p-values, logarithm of p-values, and a star rating based on the p-value ranges.
#' @param data A data frame containing the input data for analysis.
#' @param target The name of the column in the data frame that represents the group labels. The default value is "group".
#' @param feature A character vector specifying the names of the features or variables to be analyzed. If NULL, all continuous features in the data frame will be used. Default value is NULL.
#' @param feature_manipulation A logical value indicating whether feature manipulation is required. If TRUE, a custom feature manipulation function will be applied to the data. Default value is FALSE.
#'
#' @return statistical result
#' @export
#' @import dplyr
#' @import tibble
#' @author Dongqiang Zeng
#' @examples
#' # Loading TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Finding microenvironmental scores associated with Gender
#' batch_wilcoxon(data = sig_stad, target = "Gender", feature = colnames(sig_stad)[69:ncol(sig_stad)])
batch_wilcoxon<-function(data, target = "group", feature = NULL, feature_manipulation = FALSE){

  data<-as.data.frame(data)
  #change-name-of-group
  colnames(data)[which(colnames(data)==target)]<-"group"

  data<-data[!is.na(data$group), ]
  data<-data[!data$group=="", ]
  data$group<-as.character(data$group)
  group_names<- unique(data$group)

  if(is.null(feature)){
    message(">>>-- `feature` must be specified, or all continuous features will be estimated...")

    index<- menu(c("all continuous features", "selected features "), title=" >>>-- Choose features:")
    if(index == 1){
      feature<-colnames(data)
      feature<-feature[sapply(data, is.numeric)]
    }else{
      stop(">>>-- Please specify the features that you want to proceed...")
    }

  }
  feature<-feature[feature%in%colnames(data)]
  if(feature_manipulation) feature<-feature_manipulation(data = data, feature = feature, print_result = F)

  # if(!identical(group_names,c("High","Low"))) message(">>>--- `group_names` should be specified...")

  message(">>>-- Grouping information: ")
  print(table(data$group))
  ###########################################
  if(length(group_names)>2) {
    print(table(data$group))
    stop("Variable has more than two levels...")
  }

  data<-data[,c("group", feature)]
  aa<-lapply(data[,feature], function(x) wilcox.test(x ~ data[, "group"], var.equal = F))
  result_mean<-data %>% dplyr:: group_by(.$group) %>%
    dplyr:: summarise_if(is.numeric, mean)

  rownames(result_mean)<-NULL
  result_mean<-result_mean %>%
    tibble:: column_to_rownames(.,var = ".$group") %>%
    base:: as.data.frame() %>%
    t() %>%
    base:: as.data.frame() %>%
    tibble:: rownames_to_column(.,var = "sig_names")

  result_mean<-as.data.frame(result_mean)
  group_names<-group_names[order(group_names)]
  colnames(result_mean)[2:3]<-group_names
  result_mean$statistic<- result_mean[,2] - result_mean[,3]

  cc<-data.frame(sig_names=feature,
                 p.value = sapply(aa, getElement, name = "p.value"))

  cc<-cc %>%  full_join(result_mean,by = "sig_names") %>%
    dplyr:: arrange(p.value) %>%
    dplyr:: mutate(p.adj = p.adjust(.$p.value,method = "BH") ) %>%
    dplyr:: mutate(log10pvalue = log10(.$p.value)* -1) %>%
    dplyr:: mutate(stars = cut(.$p.value, breaks=c(-Inf, 0.0001, 0.001, 0.01, 0.05,0.5, Inf),
                       label=c("****","***", "**", "*", "+","")))
  cc<-tibble::as_tibble(cc)
  return(cc)

}

