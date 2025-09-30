##' Calculate siganture score using PCA function
##'
##'
##' @param eset normalized count matrix; rows are all genes in the signature that shall be summarized into one score; columns are samples
##' @param methods character vector defining whether signature scores shall be based on principal component 1 ("PC", default) or z-scores (other value)
##'
##' @return numeric vector of length ncol(eset); a score summarizing the rows of eset
##' @author Dorothee Nickles, Dongqiang Zeng
##' @export
sigScore <- function(eset, methods = "PCA") {
  if (methods == "PCA") {
    # message(paste0("Calculating siganture score using PCA function"))
    pc <- prcomp(t(eset), na.action = na.omit, scale. = T)
    sigs <- pc$x[, 1] * sign(cor(pc$x[, 1], colMeans(eset)))
  } else {
    # message(paste0("Calculating siganture score using mean of signature genes"))
    sigs <- colMeans(eset)
  }
  return(sigs)
}
#####################################
