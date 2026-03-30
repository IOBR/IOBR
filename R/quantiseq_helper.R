#' Helper Functions for quanTIseq
#'
#' @description
#' Helper functions for quanTIseq deconvolution method.
#' Source code adapted from https://github.com/FFinotello/quanTIseq.
#'
#' @name quantiseq_helper
NULL

#' Fix Expression Mixture Matrix
#'
#' @description
#' Processes expression matrix by mapping genes, converting log values,
#' quantile normalization (if arrays), and TPM normalization.
#'
#' @param mix.mat Expression matrix with genes as rows.
#' @param arrays Logical indicating if data is from arrays. Default is FALSE.
#'
#' @return Processed expression matrix.
#'
#' @keywords internal
fixMixture <- function(mix.mat, arrays = FALSE) {
  # Map gene names
  mix.mat <- mapGenes(mix.mat)

  # Un-log data in log2 base
  if (max(mix.mat) < 50) {
    mix.mat <- 2^mix.mat
  }

  # Quantile normalization
  if (arrays) mix.mat <- makeQN(mix.mat)

  # TPM normalization
  mix.mat <- t(t(mix.mat) * 1e6 / apply(mix.mat, 2, sum))

  mix.mat
}

#' Quantile Normalization
#'
#' @param mix.mat Expression matrix.
#'
#' @return Quantile normalized matrix.
#'
#' @keywords internal
makeQN <- function(mix.mat) {
  cnames <- colnames(mix.mat)
  rnames <- rownames(mix.mat)
  rlang::check_installed("preprocessCore")
  mix.mat <- preprocessCore::normalize.quantiles(as.matrix(mix.mat))
  colnames(mix.mat) <- cnames
  rownames(mix.mat) <- rnames
  mix.mat
}

#' Map Gene Names to Approved Symbols
#'
#' @description
#' Maps gene symbols to approved HGNC symbols, handling withdrawn and
#' previous symbols.
#'
#' @param mydata Expression matrix with gene symbols as rownames.
#'
#' @return Matrix with mapped gene symbols.
#'
#' @keywords internal
mapGenes <- function(mydata) {
  HGNC <- quantiseq_data$HGNC_genenames_20170418

  curgenes <- rownames(mydata)
  newgenes <- rep(NA, length(curgenes))
  newgenes2 <- rep(NA, length(curgenes))
  ind <- match(curgenes, HGNC$ApprovedSymbol)

  # Current symbols and withdrawn ones
  genes.ind.notNA <- which(!is.na(ind))
  for (i in genes.ind.notNA) {
    if (HGNC$Status[ind[i]] == "Approved") {
      newgenes[i] <- curgenes[i]
    } else if (HGNC$Status[ind[i]] == "EntryWithdrawn") {
      next
    } else {
      Wstring <- "symbolwithdrawn,see"
      newsymbol <- gsub(Wstring, "", HGNC$ApprovedName[ind[i]])
      newgenes2[i] <- newsymbol
    }
  }

  # Not found as symbols - check previous symbols and synonyms
  genes.ind.NA <- which(is.na(ind))
  for (i in genes.ind.NA) {
    genei <- curgenes[i]

    # Previous symbol?
    ind1 <- grep(genei, HGNC$PreviousSymbols)
    for (i1 in ind1) {
      array1 <- unlist(strsplit(as.character(HGNC$PreviousSymbols[i1]), ","))
      if (any(array1 == genei)) {
        newgenes2[i] <- as.character(HGNC$ApprovedSymbol[i1])
      }
    }
    # Synonym?
    ind2 <- grep(genei, HGNC$Synonyms)
    for (i2 in ind2) {
      array2 <- unlist(strsplit(as.character(HGNC$Synonyms[i2]), ","))
      if (any(array2 == genei)) {
        newgenes2[i] <- as.character(HGNC$ApprovedSymbol[i2])
      }
    }
  }

  newgenes2[which(newgenes2 %in% setdiff(newgenes, NA))] <- NA
  ind <- intersect(which(is.na(newgenes)), which(!is.na(newgenes2)))
  newgenes[ind] <- newgenes2[ind]

  mydata <- mydata[which(!is.na(newgenes)), , drop = FALSE]
  newgenes <- newgenes[which(!is.na(newgenes))]

  # Take the median if duplicates are present
  outdata <- stats::aggregate(mydata, by = list(newgenes), FUN = stats::median)
  rownames(outdata) <- outdata[, 1]
  outdata <- outdata[, -1, drop = FALSE]
  outdata <- as.data.frame(outdata)

  outdata
}

#' Run quanTIseq Deconvolution
#'
#' @param currsig Signature matrix.
#' @param currmix Mixture matrix.
#' @param scaling Scaling factors.
#' @param method Deconvolution method: "lsei", "hampel", "huber", or "bisquare".
#'
#' @return Deconvolution results matrix.
#'
#' @keywords internal
quanTIseq <- function(currsig, currmix, scaling, method) {
  method <- rlang::arg_match(method, c("lsei", "hampel", "huber", "bisquare"))

  cgenes <- intersect(rownames(currsig), rownames(currmix))
  currsig <- as.matrix(currsig[cgenes, , drop = FALSE])
  currmix <- as.matrix(currmix[cgenes, , drop = FALSE])

  if (method == "lsei") {
    # Run deconvolution with constrained least squares
    G <- diag(ncol(currsig))
    G <- rbind(G, rep(-1, ncol(G)))
    H <- c(rep(0, ncol(currsig)), -1)

    results <- apply(currmix, 2, DClsei,
      A = currsig, G = G, H = H,
      scaling = scaling
    )
  } else {
    # Run deconvolution with robust regression
    results <- apply(currmix, 2, DCrr,
      A = currsig,
      method = method, scaling = scaling
    )
  }

  t(results)
}

#' Constrained Least Squares Deconvolution
#'
#' @param b Observation vector.
#' @param A Signature matrix.
#' @param G Constraint matrix.
#' @param H Constraint vector.
#' @param scaling Scaling factor.
#'
#' @return Estimated cell fractions.
#'
#' @keywords internal
DClsei <- function(b, A, G, H, scaling) {
  rlang::check_installed("limSolve")
  sc <- norm(A, "2")
  A <- A / sc
  b <- b / sc

  res <- limSolve::lsei(
    A = A,
    B = b,
    G = G,
    H = H,
    verbose = FALSE
  )
  est <- res$X

  est.sum <- sum(est)
  est <- est / scaling
  est <- est / sum(est) * est.sum
  est <- c(est, pmax(0, 1 - sum(est)))
  names(est)[length(est)] <- "Other"

  est
}

#' Robust Regression Deconvolution
#'
#' @param b Observation vector.
#' @param A Signature matrix.
#' @param method Robust regression method.
#' @param scaling Scaling factor.
#'
#' @return Estimated cell fractions.
#'
#' @keywords internal
DCrr <- function(b, A, method, scaling) {
  rlang::check_installed("MASS")
  m <- paste0("psi.", method)

  if (m == "psi.hampel") {
    bres <- MASS::rlm(b ~ A, psi = m, a = 1.5, b = 3.5, c = 8, maxit = 1e3)
  } else {
    bres <- MASS::rlm(b ~ A, psi = m, maxit = 1e3)
  }
  est <- bres$coefficients

  # Remove intercept
  est <- est[-1]

  # Set negative values to 0 and normalize
  est[est < 0] <- 0
  est <- est / sum(est)

  # Scale by mRNA content
  est.sum <- sum(est)
  est <- est / scaling
  est <- est / sum(est) * est.sum

  names(est) <- gsub("^A", "", names(est))

  est
}
