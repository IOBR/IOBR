



#' Title
#'
#' @param eset
#' @param data_type
#' @param methods
#' @param prop
#' @param quantile
#' @param min.mad
#' @param feas
#'
#' @return
#' @export
#'
#' @examples
find_variable_genes <- function(eset, data_type = c("count", "normalized"), methods = c( "low", "mad"), prop = 0.7,
                                quantile = c(0.75, 0.5, 0.25), min.mad = 0.1, feas = NULL){


  eset <- as.matrix(eset)
  feas0 <- feature_manipulation(data = eset, is_matrix = TRUE)
  eset <- eset[feas0, ]

  if(data_type=="count" & c("low")%in%methods){
    message(paste0(">>>== Genes expressed in more than", prop*100, " percent of the samples") )
    print(table(rowSums(eset==0) < ncol(eset)*prop))
    keep<- rowSums(eset==0) < ncol(eset)*prop
    feas1 <- rownames(eset)[keep]
  }else{
    feas1 <- NULL
  }

  if("mad"%in%methods){

    eset <- log2eset(eset)
    message(paste0(">>>== min.mad = ", min.mad) )
    m.mad <- apply(eset,1,mad)
    message(paste0(">>>== Range of mad: ") )
    print(range(m.mad))

    if(quantile==0.75) index = 4
    if(quantile==0.50) index = 3
    if(quantile==0.25) index = 2
    feas2 <- rownames(eset)[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[3], min.mad))]
  }else{
    feas2 <- NULL
  }
  ####################################
  feas <- unique(c(feas1, feas2, feas))
  feas <- feas[!is.na(feas)]
  #####################################
  fx <- eset[rownames(eset)%in%feas, ]
  return(fx)
}
