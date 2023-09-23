


#' Batch to conduct correlation analysis with multiple features
#'
#' @param data A data frame containing the dataset.
#' @param target A character specifying the name of the target variable.
#' @param feature A character vector specifying the names of the feature variables.
#' @param method A character specifying the correlation method to be used. Default value is "spearman".
#'
#' @return
#' @export
#'
#' @examples
#'
#' # Loading TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Finding micro environmental scores associated with CD8 T cells
#' batch_cor(data = sig_stad, target = "CD_8_T_effector", feature = colnames(sig_stad)[69:ncol(sig_stad)])

batch_cor<-function(data, target, feature, method = "spearman"){

  data<-as.data.frame(data)
  feature<-feature[feature%in%colnames(data)]
  feature<-feature_manipulation(data = data, feature = feature)

  feature<-feature[!feature==target]

  aa<-lapply(data[,feature], function(x) cor.test(x,data[,target],method = method))

  bb<-lapply(data[,feature], function(x) exact_pvalue(x,data[,target],method = method))

  bb<-as.data.frame(bb)
  bb<-as.data.frame(t(bb))

  cc<-data.frame(sig_names=feature,
                 p.value = bb$V1,
                 statistic = sapply(aa, getElement, name = "estimate"))
  cc<-cc[order(cc$p.value, decreasing = FALSE),]
  cc$p.adj<-p.adjust(cc$p.value,method = "BH")
  cc$log10pvalue<--1*log10(cc$p.value)
  rownames(cc)<-NULL
  cc$stars <- cut(cc$p.value, breaks=c(-Inf,0.0001, 0.001, 0.01, 0.05,0.5, Inf),
                  label=c("****","***", "**", "*", "+",""))
  cc<-tibble::as_tibble(cc)
  return(cc)
}


#' Calculate exact p value of correlation
#'
#' @param x variables
#' @param y variables
#' @param method method used to conduct correlation analysis
#'
#' @return
#' @export
#' @example
#'
#' # Loading TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' exact_pvalue(sig_stad$CD8.T.cells, sig_stad$CD_8_T_effector, method = "spearman")

exact_pvalue<-function(x,y,method){

  l <- length(y)
  r <- cor(x = x, y = y, method = method, use = "complete.obs")

  if(r<0){
    t <- r * sqrt((l - 2) / (1 - r^2))
    s <- (l^3 - l) * (1 - r) / 6
    p <- 2 * (pt(q = t, df = l - 2))
  }else{
    t <- r * sqrt((l - 2) / (1 - r^2))
    s <- (l^3 - l) * (1 - r) / 6
    p <- 2 * (pt(q = -t, df = l - 2))
  }
  return(p)
}

