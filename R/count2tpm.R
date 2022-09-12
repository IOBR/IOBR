





#' Convert read counts to transcripts per million (TPM)
#'
#' @docType methods
#' @name count2tpm
#' @rdname count2tpm
#'
#' @param countMat A read count matrix, with geneid as rownames and sample as columns.
#' @param idType Type of gene id. "Ensembl" "ENTREZ","SYMBOL",
#' @param org Organism, hsa or mmu.
#' @param source default is "web", user can provide effLength manually, if idType is ensembl and source is set to "default", effLength would be provided by IOBR which was estimated by gencode.v22
#' @param effLength user can provide gene Length data with gene id and effective length
#' @param id id matched to row names of expression set
#' @param length column name of gene length
#' @param gene_symbol column name of gene symbol
#' @param remove_redundancy dafault is mean, other options are sd or median
#'
#' @return A tpm expression profile.
#'
#' @author Wubing Zhang
#' @author Dongqiang Zeng
#' @export
#'
count2tpm <- function(countMat, idType = "Ensembl", org="hsa",  source = "web", effLength = NULL, id = "id", gene_symbol = "symbol", length = "eff_length", remove_redundancy = "mean")
{
  # requireNamespace("biomaRt")
  if(!is.matrix(countMat))  countMat<-as.matrix(countMat)

  if(sum(is.na(countMat))>0){
    message(paste0("There are ", sum(is.na(countMat)) ," missing value in count matrix, these genes will be removed."))
    feas<-feature_manipulation(data = countMat, feature = rownames(countMat), is_matrix = T)
    countMat<-countMat[rownames(countMat)%in%feas,]
  }

  feas<-feature_manipulation(data = countMat, feature = rownames(countMat), is_matrix = T)
  countMat<-countMat[rownames(countMat)%in%feas,]

  if(is.null(effLength) & source == "web"){
    datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                        "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
    type = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "start_position", "end_position")
    if(org =="mmu") type[3] = "mgi_symbol"
    # listEnsemblArchives()
    # listMarts()
    # listAttributes()
    ds <- datasets[grepl(org, datasets)]
    mart <- biomaRt::useMart(host = "www.ensembl.org", biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
    ensembl <- biomaRt::getBM(attributes=type, mart = mart)
    ensembl$Length <- abs(ensembl$end_position - ensembl$start_position)

    if(toupper(idType) == "ENSEMBL"){
      len <- ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), "Length"]
      rownames(countMat) = ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), 3]
    }
    else if(toupper(idType) == "SYMBOL")
      len <- ensembl[match(rownames(countMat), ensembl[,3]), "Length"]
    else if(toupper(idType) == "ENTREZ")
      len <- ensembl[match(rownames(countMat), ensembl[,2]), "Length"]
    else
      stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")
  }


  if(source == "default" & tolower(idType) == "ensembl") {

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]

    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")

    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]

    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]

    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]

    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))

  }

  if(!is.null(effLength)){

    effLength<-as.data.frame(effLength)
    colnames(effLength)[which(colnames(effLength)==id)]<-"id"
    colnames(effLength)[which(colnames(effLength)==length)]<-"eff_length"
    effLength<-effLength[!duplicated(effLength$id),]


    countMat<-as.matrix(countMat)
    countMat<-countMat[rownames(countMat)%in%effLength$id,]
    effLength<-effLength[effLength$id%in%rownames(countMat),]

    if(id!=gene_symbol){
      # countMat<-as.matrix(countMat)
      colnames(effLength)[which(colnames(effLength)==gene_symbol)]<-"gene_symbol"
      rownames(countMat)<- effLength[match(rownames(countMat),effLength$id), "gene_symbol"]

    }else{
      # countMat<-as.matrix(countMat)
      effLength$gene_symbol<-effLength$id
      # colnames(effLength)[which(colnames(effLength)==gene_symbol)]<-"gene_symbol"
      rownames(countMat)<- effLength[match(rownames(countMat),effLength$id), "gene_symbol"]
    }

    len<- effLength[match(rownames(countMat), effLength[,"gene_symbol"]), "eff_length"]

  }

  na_idx <- which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
    countMat <- countMat[!is.na(len),]
    len = len[!is.na(len)]
  }
  tmp <- countMat / len
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))

  if(tolower(remove_redundancy)=="mean"){
    order_index <- apply(TPM,1,function(x) mean(x,na.rm=T))
  }else if(tolower(remove_redundancy)=="sd"){
    order_index <- apply(TPM,1,function(x) sd(x,na.rm=T))
  }else if(tolower(remove_redundancy)=="median"){
    order_index <- apply(TPM,1,function(x) median(x,na.rm=T))
  }

  TPM <-TPM[order(order_index,decreasing=T),]
  TPM <- TPM[!duplicated(rownames(TPM)),]
  TPM <- TPM[!is.na(rownames(TPM)),]
  TPM <- TPM[!rownames(TPM)==" ",]
  TPM <- TPM[,!is.na(colnames(TPM))]
  TPM <- TPM[,!colnames(TPM)==" "]
  return(TPM)
}
