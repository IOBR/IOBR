





#' Convert read counts to transcripts per million (TPM)
#'
#' @docType methods
#' @name count2tpm
#' @rdname count2tpm
#'
#' @param countMat A read count matrix, with geneid as rownames and sample as columns.
#' @param idType Type of gene id. "Ensembl" "ENTREZ","SYMBOL",
#' @param org Organism, hsa or mmu.
#'
#' @return A tpm expression profile.
#'
#' @author Wubing Zhang
#' @author Dongqiang Zeng
#' @export
#'
count2tpm <- function(countMat, idType = "Ensembl", org="hsa")
{
  # requireNamespace("biomaRt")
  if(class(countMat)!="matrix")  countMat<-as.matrix(countMat)
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  type = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "start_position", "end_position")
  if(org=="mmu") type[3] = "mgi_symbol"
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

  countMat <-as.matrix(countMat)
  na_idx <- which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
    countMat <- countMat[!is.na(len),]
    len = len[!is.na(len)]
  }
  tmp <- countMat / len
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))

  order_index <- apply(TPM,1,function(x) mean(x,na.rm=T))
  TPM <-TPM[order(order_index,decreasing=T),]
  TPM <- TPM[!duplicated(rownames(TPM)),]
  return(TPM)
}
