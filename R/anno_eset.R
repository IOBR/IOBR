


#' Annotating gene expression matrix and remove duplicated genes
#'
#' @description The anno_eset function is used to annotate an ExpressionSet object (eset) with gene symbols using the provided annotation data. Within the function, probes with missing symbols or symbols labelled as "NA_NA" are filtered out. The function calculates and prints the percentage of probes in the expression set that were annotated successfully. It then filters the eset to keep only the rows with probes that have matching identifiers in the annotation data. The function performs additional operations such as merging annotation data with the eset, removing unnecessary columns, transforming rows and columns, handing duplicates based on the specified method. In case of duplicates, gene symbols may be retained based on their mean values (if method is set as "mean") or standard deviation values (if method is set as "sd"). Finally, the function removes rows with all zero values, all NA values, or a NA for the first column. The output is the annotated and cleaned expression set.
#' @param eset (Required): An ExpressionSet object containing gene expression data.
#' @param annotation (Required): A data.frame that contains annotation information for the probes in the expression set.
#'  user can choose `anno_hug133plus2`, `anno_rnaseq` and `anno_illumina` as input
#' @param symbol (Optional, defaults to "symbol"): The column name in the annotation data.frame that represents the gene symbol.
#' @param probe (Optional, defaults to "probe_id"): The column name in the annotation data.frame that represents the probe identifiers.
#' @param method (Optional, defaults to "mean"): The method used to handle duplicate gene symbols; can be either "mean" or "sd".
#'
#' @return modified gene expression set
#' @export
#'
#' @examples
#' # For affymatrix data
#' data(eset_gse62254, package = "IOBR")
#' eset <- anno_eset(eset = eset_gse62254, annotation = anno_hug133plus2)
#'
#' # For RNAseq data with ensembl id
#' data(eset_stad, package = "IOBR")
#' eset <- anno_eset(eset = eset_stad, annotation = anno_grch38, probe = "id")

anno_eset<- function(eset, annotation, symbol = "symbol", probe = "probe_id", method = "mean"){

  annotation<-as.data.frame(annotation)
  colnames(annotation)[which(colnames(annotation)==symbol)]<-"symbol"
  colnames(annotation)[which(colnames(annotation)==probe)]<-"probe_id"
  annotation<-annotation[,c("probe_id","symbol")]
  annotation<-annotation[!annotation$symbol=="NA_NA",]
  annotation<-annotation[!is.na(annotation$symbol),]

  message(paste0("Row number of original eset: "))
  message(paste0(">>>>  ",dim(eset)[1]))

  annotation<-annotation[!is.na(annotation$symbol),]

  anno_count<-length(rownames(eset)[rownames(eset)%in%annotation$probe_id])/length(rownames(eset))
  message(paste0(paste0(sprintf(">>> %1.2f%%", 100*anno_count)," of probe in expression set was annotated")))

  annotation<-annotation[annotation$probe_id%in%rownames(eset),]
  eset<-eset[rownames(eset)%in%annotation$probe_id,]

  eset<-as.data.frame(eset)
  eset<-tibble:: rownames_to_column(eset,var = "id")
  eset<-merge(annotation,eset,by.x="probe_id",by.y="id",all = F)
  eset<-eset[,-1]

  eset<-as.data.frame(eset)
  rownames(eset)<-NULL

  symbol<-"symbol"
  dups <- dim(eset)[1] - length(unique(eset[,symbol]))

  if(dups==0){
    eset<-tibble:: column_to_rownames(eset,var = symbol)
  }else{
    if(method=="mean"){
      order_index=apply(eset[,setdiff(colnames(eset),symbol)],1,function(x) mean(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>%dplyr:: distinct(!!sym(symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = symbol)
    }else if(method == "sd"){
      order_index = apply(eset[,setdiff(colnames(eset),symbol)],1,function(x) sd(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>% distinct(!!sym(symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = symbol)
    }else if(method == "sum"){
      order_index = apply(eset[,setdiff(colnames(eset),symbol)],1,function(x) sum(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>% distinct(!!sym(symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = symbol)
    }
  }

  quit<-rowSums(eset==0)==ncol(eset)
  eset=eset[!quit,]

  quit<-rowSums(is.na(eset))==ncol(eset)
  eset=eset[!quit,]

  quit<-is.na(eset[,1])
  eset=eset[!quit,]

  message(paste0("Row number after filtering duplicated gene symbol: "))
  message(paste0(">>>>  ",dim(eset)[1]))

  return(eset)
}
