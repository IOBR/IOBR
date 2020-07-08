





#' Estimate the fraction of cell types using defined reference genes
#'
#' @param eset expression data with matched gene id of reference
#' @param reference immune cell gene matrix; eg lm22, lm6 or can be generate using generateRef/generateRef_rnaseq
#' @param scale_reference a logical value indicating whether the reference be scaled or not. If TRUE, the value in reference file will be centered and scaled in row direction.
#' @param arrays a logical value. If TRUE, the columns of the input data will be normalized to have the same quantiles.
#' @param method method must be set to either 'svm' or 'lsei'
#' @param perm Set permutations for statistical analysis (≥100 permutations recommended).
#' @return
#' @export
#' @import limSolve
#' @import limma
#' @import preprocessCore
#' @import e1071
#' @import parallel
#' @examples
#'
#' data("eset_crc")
#' data("lm22")
#' ## eset_crc has already been normalized, so here normalize = FALSE
#' til_svm <- deconvo_constru_sig(eset = eset_crc[, 1:50], reference = lm22, scale_reference = T, method = "svm", arrays = F)
#' til_lsei <- deconvo_constru_sig(eset= eset_crc[, 1:50], reference = lm22_ref, scale_reference = T, method = "lsei",arrays = F)
deconvo_constru_sig <- function(eset, reference, scale_reference = TRUE, method = "svm", arrays = TRUE, perm = 500){

  if (length(intersect(rownames(eset), rownames(reference))) == 0){
    stop("None identical gene between eset and reference had been found.
         Check your eset using: intersect(rownames(eset), rownames(reference))")
  }
  eset <- as.matrix(eset)
  reference <- as.matrix(reference)
  # quantile normalized
  if (arrays){
    roweset <- rownames(eset)
    coleset <- colnames(eset)
    eset <- normalize.quantiles(eset)
    rownames(eset) <- roweset
    colnames(eset) <- coleset
  }
  # scale_reference
  if (scale_reference){
    reference <- (reference - mean(reference))/sd(as.vector(reference))
  }
  Ymedian <- max(median(eset),1)
  # common eset
  common <- intersect(rownames(eset), rownames(reference))
  eset <- eset[match(common, rownames(eset)), ]
  reference <- reference[match(common, rownames(reference)), ]
  # deconvolution
  output <- matrix()
  if (method == "svm"){

    # message(paste0("\n", ">>> Running ", "cell estimation in SVM mode"))
    nulldist <- sort(doPerm(perm = perm, X = reference, Y = eset,absolute = FALSE, abs_method = 'sig.score')$dist)
    # absolute <- FALSE
    # abs_method <- 'sig.score'

    itor <- 1
    samples <- ncol(eset)
    while (itor <= samples) {
      y <- eset[, itor]
      y <- (y - mean(y))/sd(y)
      # svr
      result <- CoreAlg(X = reference, y = y, absolute = FALSE, abs_method = 'sig.score')
      w <- result$w
      mix_r <- result$mix_r
      mix_rmse <- result$mix_rmse
      if(absolute && abs_method == 'sig.score') {
        w <- w * median(Y[,itor]) / Ymedian
      }
      pvalue <-  1- which.min(abs(nulldist - mix_r))/length(nulldist)
      out <- c(w, pvalue, mix_r, mix_rmse)
      if (absolute){
        out <- c(out, sum(w))
      }
      if (itor == 1){output = out}else{output = rbind(output, out)}
      itor <- itor + 1
    }
    header <- c(colnames(reference),"P-value","Correlation","RMSE")
    if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
    colnames(output) <- header
    rownames(output) <- colnames(eset); output <- output[, -1]
    infiltration <- output
  }else if(method == 'lsei'){

    # message(paste0("\n", ">>> Running ", "cell estimation in lsei mode"))
    Numofx <- ncol(reference)
    AA <- reference
    EE <- rep(1, Numofx)
    FF <- 1
    GG <- diag(nrow=Numofx)
    HH <- rep(0, Numofx)

    out.all <- c()
    itor <- 1
    samples <- ncol(eset)
    while (itor <= samples){
      BB <- eset[, itor]
      BB <- (BB - mean(BB))/sd(BB)
      out <- lsei(AA, BB, EE, FF, GG, HH)
      out.all <- rbind(out.all, out$X)
      itor <- itor + 1
    }
    rownames(out.all) <- colnames(eset)
    infiltration <- out.all
  }else{
    stop("The method must be set to either 'svm' or 'lsei' ")
  }
  return(infiltration)
  }

