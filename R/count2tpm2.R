



#' Convert read counts to transcripts per million (TPM)
#'
#' @docType methods
#' @name count2tpm
#' @rdname count2tpm
#'
#' @param countMat A read count matrix, with `geneid` as row names and samples as columns.
#' @param idType Type of gene id. "Ensembl" "ENTREZ","SYMBOL"
#' @param org Organism, `hsa` or `mmus`
#' @param source default is "local". Other option is `biomart`. user can also provide `effLength` manually, if `idType` is `ensembl`, and source is set to `local`, `effLength` was provided by IOBR which was estimated by function `getGeneLengthAndGCContent` of EDASeq package at 2023-02-10.
#' @param effLength user can provide gene Length data with gene id and effective length
#' @param id id matched to row names of expression set
#' @param length column name of gene length
#' @param gene_symbol column name of gene symbol
#' @param remove_redundancy default is mean, other options are `sd`, `mean` or `median`
#'
#' @return A tpm expression profile.
#'
#' @author Wubing Zhang
#' @author Dongqiang Zeng
#' @export
#'
count2tpm <- function(countMat, idType = "Ensembl", org = "hsa",  source = "local", effLength = NULL, id = "id", gene_symbol = "symbol", length = "eff_length", remove_redundancy = "mean"){


  # requireNamespace("biomaRt")
  if(!is.matrix(countMat))  countMat<-as.matrix(countMat)

  if(sum(is.na(countMat))>0){
    message(paste0("There are ", sum(is.na(countMat)) ," missing value in count matrix, these genes will be removed."))
    feas<-feature_manipulation(data = countMat, feature = rownames(countMat), is_matrix = T)
    countMat<-countMat[rownames(countMat)%in%feas,]
  }

  feas<-feature_manipulation(data = countMat, feature = rownames(countMat), is_matrix = T)
  countMat<-countMat[rownames(countMat)%in%feas,]

  if(is.null(effLength) & source == "biomart"){
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
    #######################################

    ensembl$Length <- abs(ensembl$end_position - ensembl$start_position)

    message(">>>--- This function is being optimised and we strongly recommend that you should set `source` as `local`....")
    #######################################
    if(toupper(idType) == "ENSEMBL"){
      len <- ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), "Length"]
      rownames(countMat) = ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), 3]
    }
    else if(toupper(idType) == "SYMBOL")
      len <- ensembl[match(rownames(countMat), ensembl[,3]), "Length"]
    else if(toupper(idType) == "ENTREZ")
      len <- ensembl[match(rownames(countMat), ensembl[,2]), "Length"]
    else
      stop("Please input right type of gene name, such as `ensembl`, `entrez`, or `symbol` ...")
  }


  if(source == "local" & tolower(idType) == "ensembl" & org == "hsa") {

    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_grch38[,c("id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }else if(source == "local" & tolower(idType) == "entrez"  & org == "hsa"){


    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL as row names")
    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_grch38[,c("entrez", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]
    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat), ]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }else if(source == "local" & tolower(idType) == "symbol"  & org == "hsa"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL as row names")
    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")
    length_ensembl<-anno_grch38[,c("symbol", "eff_length", "gc")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]

    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]

    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 1]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }

  #######################################################################
  if(source == "local" & tolower(idType) == "ensembl" & org == "mmus") {

    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_gc_vm32[,c("id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }else if(source == "local" & tolower(idType) == "mgi"  & org == "mmus"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL ID as row names")
    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_gc_vm32[,c("mgi_id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]
    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat) <- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }else if(source == "local" & tolower(idType) == "symbol"  & org == "mmus"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL ID as row names")
    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")
    length_ensembl<-anno_gc_vm32[,c("symbol", "eff_length", "gc")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 1]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }

  #########################################################################
  if(!is.null(effLength)){
    effLength<-as.data.frame(effLength)
    colnames(effLength)[which(colnames(effLength)==id)]<-"id"
    colnames(effLength)[which(colnames(effLength)==length)]<-"eff_length"
    effLength<-effLength[!duplicated(effLength$id),]

    countMat<-as.matrix(countMat)
    countMat<-countMat[rownames(countMat)%in%effLength$id, ]
    effLength<-effLength[effLength$id%in%rownames(countMat), ]

    if(id!= gene_symbol){
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
    warning(paste0(">>>--- Omit ", length(na_idx), " genes of which length is not available !"))
    countMat <- countMat[!is.na(len), ]
    len = len[!is.na(len)]
  }
  #####################################
  tmp <- countMat / c(len/1000) # (`per million` scaling factor)
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
  # print(head(TPM))
  # rownames(TPM) <- as.character(rownames(TPM))
  # TPM <- BiocGenerics:: as.data.frame(TPM)
  # print(head(TPM))
  #####################################
  # TPM <- rownames_to_column(TPM, var = "symbol" )
  # TPM <- remove_duplicate_genes(eset = TPM, column_of_symbol = "symbol", method = tolower(remove_redundancy))

  # if(tolower(remove_redundancy)=="mean"){
  #   order_index <- apply(TPM,1,function(x) mean(x,na.rm=T))
  # }else if(tolower(remove_redundancy)=="sd"){
  #   order_index <- apply(TPM,1,function(x) sd(x,na.rm=T))
  # }else if(tolower(remove_redundancy)=="median"){
  #   order_index <- apply(TPM,1,function(x) median(x,na.rm=T))
  # }
  # TPM <- TPM[order(order_index,decreasing=T),]

  # TPM <- TPM[!duplicated(rownames(TPM)),]
  TPM <- TPM[!is.na(rownames(TPM)),]
  TPM <- TPM[!rownames(TPM)==" ",]

  TPM <- rownames_to_column(as.data.frame(TPM), var = "symbol")
  TPM <- remove_duplicate_genes(eset = TPM, column_of_symbol = "symbol")
  # TPM <- TPM[,!is.na(colnames(TPM))]
  # TPM <- TPM[,!colnames(TPM)==" "]
  return(TPM)
}


# help("count2tpm")
# eset<-load_rna(project = "TCGA", ProjectID = "TCGA-STAD", tissue_type = "tumor_normal", data_type = "count")
# mhead(eset)
# rownames(eset)<-substring(rownames(eset), 1, 15)
# eset_tpm<-count2tpm( countMat = eset,   idType = "Ensembl",   org = "hsa",  source = "web")
######################################
# get gene data from BioMart
# ensembl<- EDASeq:: getGeneLengthAndGCContent(id = anno_grch38$ensgene, org = "hsa", mode = "biomart")
# ensembl<- rownames_to_column(ensembl, var = "id")
# ensembl<- merge(ensembl, anno_grch38, by.x = "id", by.y = "ensgene", all = FALSE)
#######################################
# https://cloud.tencent.com/developer/article/2031969




