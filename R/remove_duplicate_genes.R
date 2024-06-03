


#' Remove duplicate gene symbol on a gene expression data
#'
#' This function removes duplicate gene symbols from a gene expression dataset. It offers the option to aggregate duplicate gene expressions using different methods such as mean or standard deviation.
#'
#' @param eset Gene expression data with gene symbol in `column_of_symbol`
#' @param column_of_symbol name of column contained gene symbols
#' @param method method used to filter duplicate genes; default is mean value
#' @import tibble dplyr
#' @importFrom rlang sym
#' 
#' @return A modified gene expression data frame with duplicates removed and selected method applied for aggregation.
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#'
#' # loading eset
#' data("eset_stad", package = "IOBR")
#' # annotation
#' eset_stad <- anno_eset(eset = eset_stad, annotation = anno_rnaseq)
#' eset_stad <- rownames_to_column(eset_stad, var = "id")
#'
#' # Creating duplicate gene names
#' eset_stad[2:3, "id"] <- "MT-CO1"
#' # Counting the number of identical names
#' summary(duplicated(eset_stad$id))
#' # De-duplication of rows with the same gene name using the average value
#' eset_stad<-remove_duplicate_genes(eset = eset_stad, column_of_symbol = "id", method = "mean")
#' summary(duplicated(eset_stad$id))

remove_duplicate_genes<-function(eset, column_of_symbol, method = "mean"){
  eset<-as.data.frame(eset)
  rownames(eset)<-NULL

  dups <- dim(eset)[1] - length(unique(eset[,column_of_symbol]))

  if(dups==0){
    eset<-tibble:: column_to_rownames(eset,var = column_of_symbol)
    return(eset)
  }else{
    if(method=="mean"){
      order_index=apply(eset[,setdiff(colnames(eset),column_of_symbol)],1,function(x) mean(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>%dplyr:: distinct(!!sym(column_of_symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = column_of_symbol)
      return(eset)
    }else if(method == "sd"){
      order_index = apply(eset[,setdiff(colnames(eset),column_of_symbol)],1,function(x) sd(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>% distinct(!!sym(column_of_symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = column_of_symbol)
      return(eset)
    }
  }
}

