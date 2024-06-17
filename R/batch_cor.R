


#' Batch Correlation Analysis
#'
#' Performs a batch correlation analysis between a target variable and multiple feature variables in a dataset.
#' This function allows for the selection of the correlation method and applies corrections for multiple testing.
#'
#' @param data A data frame containing the dataset with both target and feature variables.
#' @param target A character string specifying the name of the target variable within the dataset.
#' @param feature A character vector specifying the names of the feature variables to be analyzed.
#' @param method A character string specifying the correlation method to be used; default is "spearman".
#'
#' @return Returns a tibble containing the feature names, p-values, correlation coefficients, adjusted p-values,
#'         log-transformed p-values, and significance stars based on p-value thresholds.
#'
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Conduct correlation analysis
#' results <- batch_cor(data = sig_stad, target = "CD_8_T_effector", feature = colnames(sig_stad)[69:ncol(sig_stad)])
batch_cor<-function(data, target, feature, method = "spearman"){

  if(!target%in%colnames(data)) stop(">>>== target was not in the colnames of data")
  data<-as.data.frame(data)
  feature<-feature[feature%in%colnames(data)]

  if(length(feature)==0) stop(">>>== features were not in the colnames of data")

  feature<-feature_manipulation(data = data, feature = feature)

  feature<-feature[!feature==target]


  aa<-lapply(data[,feature], function(x) cor.test(x,data[,target],method = method))

  bb<-lapply(data[,feature], function(x) exact_pvalue(x,data[,target],method = method))

  bb<-as.data.frame(bb)

  bb<-as.data.frame(t(bb))
  # print(head(bb))
  cc<-data.frame(sig_names = feature,
                 p.value = bb$V1,
                 statistic = sapply(aa, getElement, name = "estimate"))

  # cc<-cc[base::order(cc$p.value, decreasing = FALSE),]
  cc$p.adj <- p.adjust(cc$p.value,method = "BH")
  cc$log10pvalue<- -1*log10(cc$p.value)
  rownames(cc)<-NULL
  cc$stars <- cut(cc$p.value, breaks=c(-Inf,0.0001, 0.001, 0.01, 0.05,0.5, Inf),
                  label=c("****","***", "**", "*", "+",""))
  cc<-cc[base::order(cc$p.value, decreasing = FALSE),]
  cc<-tibble::as_tibble(cc)
  print(head(cc))
  return(cc)
}


#' Calculate Exact P-value for Correlation
#'
#' Computes the exact p-value for the correlation between two variables based on a specified method.
#' This function is typically used to support detailed statistical tests in correlation studies.
#'
#' @param x A numeric vector of data points corresponding to the first variable.
#' @param y A numeric vector of data points corresponding to the second variable.
#' @param method A character string specifying the correlation method to be used; supports "spearman", "pearson", etc.
#'
#' @return Returns a single numeric value representing the p-value for the correlation test.
#'
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Calculate exact p-value for correlation between two variables
#' p_val <- exact_pvalue(x = sig_stad$CD8.T.cells, y = sig_stad$CD_8_T_effector, method = "spearman")
#' print(p_val)
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

