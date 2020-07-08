


#' Annotating gene expression matrix and remove duplicated genes
#'
#' @param eset gene expression matrix with probe or ensembl, entrez id in the rownames
#' @param annotation annotation file with gene symbol and other ids matched to rowname of eset
#'  user can choose anno_hug133plus2, anno_rnaseq and anno_illumina as input
#' @param symbol column name of gene symbol in annotation file
#' @param probe column name of relevant id in annotation file
#' @param method method used to remove duplicated genes
#'
#' @return modified gene expression set
#' @export
#'
#' @examples
anno_eset<- function(eset,annotation,symbol = "symbol", probe = "probe_id",method = "mean"){

  annotation<-as.data.frame(annotation)
  colnames(annotation)[which(colnames(annotation)==symbol)]<-"symbol"
  colnames(annotation)[which(colnames(annotation)==probe)]<-"probe_id"
  annotation<-annotation[,c("probe_id","symbol")]
  annotation<-annotation[!annotation$symbol=="NA_NA",]

  message(paste0("Row number of original eset: "))
  message(paste0(">>>>  ",dim(eset)[1]))

  annotation<-annotation[!is.na(annotation$symbol),]
  annotation<-annotation[annotation$probe_id%in%rownames(eset),]
  eset<-eset[rownames(eset)%in%annotation$probe_id,]

  eset<-as.data.frame(eset)
  eset<-tibble:: rownames_to_column(eset,var = "id")
  eset<-merge(annotation,eset,by.x="probe_id",by.y="id",all = F)
  eset<-eset[,-1]

  eset<-as.data.frame(eset)
  rownames(eset)<-NULL

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
    }
  }

  quit<-rowSums(eset==0)==ncol(eset)
  eset=eset[!quit,]

  quit<-rowSums(is.na(eset))==ncol(eset)
  eset=eset[!quit,]

  quit<-is.na(eset[,1])
  eset=eset[!quit,]

  message(paste0("Row number of filtered eset: "))
  message(paste0(">>>>  ",dim(eset)[1]))

  return(eset)
}
