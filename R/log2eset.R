



#' Judge whether the expression set needs log2 transformation
#'
#' @param eset gene expression set with row as genes and sample in the column
#'
#' @return gene expression set
#' @export
#'
#' @examples
#' # Loading TCGA-STAD expresion data(raw count matrix)
#' data("eset_stad", package = "IOBR")
#' # transform count data to tpm
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' eset <- log2eset(eset)
log2eset<-function(eset = NULL){

  qx <- as.numeric(quantile(eset, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
  log_judge <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

  if(log_judge){
    eset[eset < 0] <- 0
    eset <- log2(eset+1)
    message(">>> log2 transformation was finished")
  }else{
    message(">>> log2 transformation is not necessary")}

  return(eset)
}
