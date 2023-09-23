


#' Title find_variable_genes
#'
#' @description The "find_variable_genes" function is designed to identify and filter variable genes in a gene expression matrix based on the specified criteria. It takes a gene expression matrix ("x") as input and returns a filtered matrix ("fx") containing only the variable genes.
#' @param x A gene expression matrix. It should be a numeric matrix or a matrix-like object.
#' @param data_type The type of data in the gene expression matrix. The default value is "count", indicating count data. This parameter determines the method used for filtering variable genes.
#' @param prop The proportion threshold for low expression filtering. It specifies the minimum proportion of samples in which a gene should be expressed to be considered as a variable gene. The default value is 0.7, meaning genes expressed in more than 7% of the samples will be kept.
#' @param methods A character vector specifying the filtering methods to be applied. The default value is c("low", "mad"), indicating both low expression filtering and median absolute deviation (MAD) filtering methods.
#' @param quantile The quantile value for MAD filtering. It determines the threshold for selecting genes based on their median absolute deviation. The default value is 0.75, indicating the top 75% of genes with higher MAD will be selected.
#' @param min.mad The minimum median absolute deviation threshold for MAD filtering. It specifies the minimum MAD value a gene should have to be considered as a variable gene. The default value is 0.1.
#' @param feas An optional vector of gene names or indices specifying additional genes to be considered as variable genes. This allows users to manually include certain genes in the filtered matrix. The default value is NULL.
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#' @examples
find_variable_genes <- function(x, data_type = c("count", "normalized"), methods = c( "low", "mad"), prop = 0.7, quantile = c(0.75, 0.5, 0.25), min.mad = 0.1, feas = NULL){


  x <- as.matrix(x)
  feas0 <- feature_manipulation(data = x, is_matrix = TRUE)
  x <- x[feas0, ]

  if(data_type=="count" & c("low")%in%methods){
    # 去掉低表达的基因

    message(paste0(">>>== Genes expressed in more than", prop*100, "% of the samples") )
    print(table(rowSums(x==0) < ncol(x)*prop))
    keep<- rowSums(x==0) < ncol(x)*prop
    feas1 <- rownames(x)[keep]
  }else{
    feas1 <- NULL
  }

  if("mad"%in%methods){

    x <- log2eset(x)
    message(paste0(">>>== min.mad = ", min.mad) )
    m.mad <- apply(x,1,mad)
    message(paste0(">>>== Range of mad: \n") )
    print(range(m.mad))

    if(quantile==0.75) index = 4
    if(quantile==0.50) index = 3
    if(quantile==0.25) index = 2
    feas2 <- rownames(x)[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[3], min.mad))]
  }else{
    feas2 <- NULL
  }
  ####################################
  feas <- unique(c(feas1, feas2, feas))
  feas <- feas[!is.na(feas)]
  #####################################
  fx <- x[rownames(x)%in%feas, ]
  return(fx)
}
