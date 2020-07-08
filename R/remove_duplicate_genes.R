


#' Remove duplicate gene symbol on a gene expression data
#'
#' @param eset Gene expression data with gene symbol in `column_of_symbol`
#' @param column_of_symbol name of column contained gene symbols
#' @param method method used to filter duplicate genes; default is mean value
#' @import tibble
#' @return
#' @export
#'
#' @examples
#' eset[1:5,1:5]
#' eset<-remove_duplicate_genes(eset = eset,column_of_symbol = "symbol",method = "mean")
#' summary(duplicated(rownames(eset)))
#

remove_duplicate_genes<-function(eset, column_of_symbol ,method = "mean"){
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

