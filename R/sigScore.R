#' Calculate Signature Score Using PCA or Mean Methods
#'
#' Computes signature scores from gene expression data using either Principal Component Analysis (PCA) or mean-based approaches.
#'
#' @param eset Normalized count matrix; rows are genes in the signature to be summarized into one score; columns are samples.
#' @param methods Character vector defining whether signature scores shall be based on principal component 1 ("PCA", default) or mean values (any other value).
#'
#' @return Numeric vector of length ncol(eset); a score summarizing the rows of eset.
#' @author Dorothee Nickles, Dongqiang Zeng
#' @export
#' @examples
#' # Example usage with PCA method
#' sigScore(eset = expression_matrix, methods = "PCA")
#' # Example usage with mean method
#' sigScore(eset = expression_matrix, methods = "mean")
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
