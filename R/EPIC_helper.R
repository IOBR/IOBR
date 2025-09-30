# ####################################'
#
# Package EPIC - developed by Julien Racle from David Gfeller's group of the
# University of Lausanne and Ludwig Institute for Cancer Research. Copyrights
# remain reserved as described in the included LICENSE file.
#
# ####################################'

#' Estimate the proportion of immune and cancer cells.
#'
#' \code{EPIC} takes as input bulk gene expression data (RNA-seq) and returns
#'  the proportion of mRNA and cells composing the various samples.
#'
#' This function uses a constrained least square minimization to estimate the
#' proportion of each cell type with a reference profile and another
#' uncharacterized cell type in bulk gene expression samples.
#'
#' The names of the genes in the bulk samples, the reference samples and in the
#' gene signature list need to be the same format (gene symbols are used in the
#' predefined reference profiles). The full list of gene names don't need to be
#' exactly the same between the reference and bulk samples: \emph{EPIC}
#' will use the intersection of the genes. In case of duplicate gene names,
#' \emph{EPIC} will use the median value per duplicate - if you want to consider
#' these cases differently, you can remove the duplicates before calling
#' \emph{EPIC}.
#'
#' @param bulk A matrix (\code{nGenes} x \code{nSamples}) of the genes
#'    expression from each bulk sample (the counts should be given in TPM,
#'    RPKM or FPKM when using the prebuilt reference profiles). This matrix
#'    needs to have rownames telling the gene names (corresponds to the gene
#'    symbol in the prebuilt reference profiles (e.g. CD8A, MS4A1) - no
#'    conversion of IDs is performed at the moment).
#'    It is advised to keep all genes in the bulk instead of a subset of
#'    signature genes (except if \code{scaleExprs = FALSE} in which case it
#'    doesn't make any difference).
#' @param reference (optional): A string or a list defining the reference cells.
#'    It can take multiple formats, either: \itemize{
#'      \item\code{NULL}: to use the default reference profiles and genes
#'        signature \code{\link{TRef}}.
#'      \item a char: either \emph{"BRef"} or \emph{"TRef"}
#'        to use the reference cells and genes signature of the corresponding
#'        datasets (see \code{\link{BRef}} and \code{\link{TRef}}).
#'      \item a list. When a list it should include: \describe{
#'        \item{\code{$refProfiles}}{a matrix (\code{nGenes} x \code{nCellTypes})
#'        of the reference cells genes expression (don't include a column of
#'        the 'other cells' (representing usually the cancer cells for which
#'        such a profile is usually not conserved between samples);
#'        the rownames needs to be defined as well as the colnames giving the
#'        names of each gene and reference cell types respectively.
#'        It is advised to keep all genes in this \code{refProfiles} matrix
#'        instead of a subset of signature genes;
#'        }
#'        \item{\code{$sigGenes}}{a character vector of the gene names to use as
#'          signature - sigGenes can also be given as a direct input to EPIC
#'          function;}
#'        \item{\code{$refProfiles.var}}{(optional): a matrix (\code{nGenes} x
#'        \code{nCellTypes}) of the variability of each gene expression for each
#'        cell type, which is used to define weights on each gene for the
#'        optimization (if this is absent, we assume an identical variability
#'        for all genes in all cells) - it needs to have the same dimnames than
#'        refProfiles;}
#'        }
#'    }
#' @param mRNA_cell (optional): A named numeric vector: tells (in arbitrary
#'    units) the amount of mRNA for each of the reference cells and of the
#'    other uncharacterized (cancer) cell. Two names are of special meaning:
#'    \emph{"otherCells"} - used for the mRNA/cell value of the "other cells"
#'    from the sample (i.e. the cell type that don't have any reference gene
#'    expression profile) ; and \emph{default} - used for the mRNA/cell of the
#'    cells from the reference profiles for which no specific value is given
#'    in mRNA_cell (i.e. if mRNA_cell=c(Bcells=2, NKcells=2.1, otherCells=3.5,
#'    default=1), then if the refProfiles described Bcells, NKcells and Tcells,
#'    we would use a value of 3.5 for the "otherCells" that didn't have any
#'    reference profile and a default value of 1 for the Tcells when computing
#'    the cell fractions).
#'    To note: if data is in tpm, this mRNA per cell would ideally correspond
#'    to some number of transcripts per cell.
#' @param mRNA_cell_sub (optional): This can be given instead of \code{mRNA_cell} (or
#'    in addition to it). It is also a named numeric vector, used to replace
#'    only the mRNA/cell values from some cell types (or to add values for new
#'    cell types). The values given in mRNA_cell_sub will overwrite the default
#'    values as well as those that might have been given by mRNA_cell.
#' @param sigGenes (optional): a character vector of the gene names to use as
#'    signature for the deconvolution. In principle this is given with the
#'    reference as the "reference$sigGenes" but if we give a value for this
#'    input variable, it is these signature genes that will be used instead of
#'    the ones given with the reference profile.
#' @param scaleExprs (optional, default is TRUE): boolean telling if the bulk
#'    samples and reference gene expression profiles should be rescaled based on
#'    the list of genes in common between the them (such a rescaling is
#'    recommanded).
#' @param withOtherCells (optional, default is TRUE): if EPIC should allow for
#'    an additional cell type for which no gene expression reference profile is
#'    available or if the bulk is assumed to be composed only of the cells with
#'    reference profiles.
#' @param constrainedSum (optional, default is TRUE): tells if the sum of all
#'    cell types should be constrained to be < 1. When
#'    \code{withOtherCells=FALSE}, there is additionally a constrain the the sum
#'    of all cell types with reference profiles must be > 0.99.
#' @param rangeBasedOptim (optional): when this is FALSE (the default), the
#'    least square optimization is performed as described in
#'    \href{https://elifesciences.org/articles/26476}{
#'     \cite{Racle et al., 2017, eLife}}, which is recommanded.
#'    When this variable is TRUE, EPIC uses the variability of each gene
#'    from the reference profiles in another way: instead of defining weights
#'    (based on the variability) for the fit of each gene, we define a range of
#'    values accessible for each gene (based on the gene expression value in
#'    the reference profile +/- the variability values). The
#'    error that the optimization tries to minimize is by how much
#'    the predicted gene expression is outside of this allowed range of values.
#'
#' @return A list of 3 matrices:\describe{
#'  \item{\code{mRNAProportions}}{(\code{nSamples} x (\code{nCellTypes+1})) the
#'    proportion of mRNA coming from all cell types with a ref profile + the
#'    uncharacterized other cell.}
#'  \item{\code{cellFractions}}{(\code{nSamples} x (\code{nCellTypes+1})) this
#'    gives the proportion of cells from each cell type after accounting for
#'    the mRNA / cell value.}
#'  \item{\code{fit.gof}}{(\code{nSamples} x 12) a matrix telling the quality
#'    for the fit of the signature genes in each sample. It tells if the
#'    minimization converged, and other info about this fit comparing the
#'    measured gene expression in the sigGenes vs predicted gene expression in
#'    the sigGenes.}
#' }
#'
#' @examples
#' res1 <- EPIC(melanoma_data$counts)
#' res1$cellFractions
#' res2 <- EPIC(melanoma_data$counts, TRef)
#' res3 <- EPIC(bulk = melanoma_data$counts, reference = TRef)
#' res4 <- EPIC(melanoma_data$counts, reference = "TRef")
#' res5 <- EPIC(melanoma_data$counts, mRNA_cell_sub = c(Bcells = 1, otherCells = 5))
#' # Various possible ways of calling EPIC function. res 1 to 4 should
#' # give exactly the same outputs, and the elements res1$cellFractions
#' # should be equal to the example predictions found in
#' # melanoma_data$cellFractions.pred for these first 4 results.
#' # The values of cellFraction for res5 will be different due to the use of
#' # other mRNA per cell values for the B and other cells.
#'
#' @export
EPIC <- function(bulk, reference = NULL, mRNA_cell = NULL, mRNA_cell_sub = NULL,
                 sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE,
                 constrainedSum = TRUE, rangeBasedOptim = FALSE) {
  # Checking the correct format of the bulk sample input
  if (!is.matrix(bulk) && !is.data.frame(bulk)) {
    stop("'bulk' needs to be given as a matrix or data.frame")
  }

  # First get the value of the reference profiles depending on the input
  # 'reference'.
  with_w <- TRUE
  if (is.null(reference)) {
    reference <- IOBR::TRef
  } else if (is.character(reference)) {
    if (reference %in% c("TRef", "BRef")) {
      reference <- get(reference, pos = "package:IOBR")
      # Replace the char defining the reference name by the corresponding
      # pre-built reference values.
    } else {
      stop(
        "The reference, '", reference, "' is not part of the allowed ",
        "references:", paste(c("TRef", "BRef"), collapse = ", ")
      )
    }
  } else if (is.list(reference)) {
    refListNames <- names(reference)
    if ((!all(c("refProfiles", "sigGenes") %in% refListNames)) ||
      (("refProfiles" %in% refListNames) && !is.null(sigGenes))) {
      stop(
        "Reference, when given as a list needs to contain at least the ",
        "fields 'refProfiles' and 'sigGenes' (sigGenes could also be ",
        "given as input to EPIC instead)"
      )
    }
    if (!is.matrix(reference$refProfiles) && !is.data.frame(reference$refProfiles)) {
      stop("'reference$refProfiles' needs to be given as a matrix or data.frame")
    }
    if (!("refProfiles.var" %in% refListNames)) {
      warning(
        "'refProfiles.var' not defined; using identical weights ",
        "for all genes"
      )
      with_w <- FALSE
    } else if (!is.matrix(reference$refProfiles.var) &&
      !is.data.frame(reference$refProfiles.var)) {
      stop(
        "'reference$refProfiles.var' needs to be given as a matrix or ",
        "data.frame when present."
      )
    } else if (!identical(dim(reference$refProfiles.var), dim(reference$refProfiles)) ||
      !identical(
        dimnames(reference$refProfiles.var),
        dimnames(reference$refProfiles)
      )) {
      stop(
        "The dimensions and dimnames of 'reference$refProfiles' and ",
        "'reference$refProfiles.var' need to be the same"
      )
    }
  } else {
    stop("Unknown format for 'reference'")
  }

  bulk <- merge_duplicates(bulk, in_type = "bulk samples")
  refProfiles <- merge_duplicates(reference$refProfiles,
    in_type = "reference profiles"
  )
  if (with_w) {
    refProfiles.var <- merge_duplicates(reference$refProfiles.var, warn = F)
    # Don't warn here as we're already warning for refProfile and they had same
    # dim names.
  } else {
    refProfiles.var <- 0
  }

  nSamples <- NCOL(bulk)
  samplesNames <- colnames(bulk)
  if (is.null(samplesNames)) {
    samplesNames <- 1:nSamples
    colnames(bulk) <- samplesNames
  }
  nRefCells <- NCOL(refProfiles)
  refCellsNames <- colnames(refProfiles)


  # Keeping only common genes and normalizing the counts based on these common
  # genes
  bulk_NA <- apply(is.na(bulk), MARGIN = 1, FUN = all)
  if (any(bulk_NA)) {
    warning(sum(bulk_NA), " genes are NA in all bulk samples, removing these.")
    bulk <- bulk[!bulk_NA, ]
  }
  bulkGenes <- rownames(bulk)
  refGenes <- rownames(refProfiles)
  commonGenes <- intersect(bulkGenes, refGenes)

  if (is.null(sigGenes)) {
    sigGenes <- unique(reference$sigGenes)
  } # Keep only once in case of duplicates
  sigGenes <- sigGenes[sigGenes %in% commonGenes]
  nSigGenes <- length(sigGenes)
  if (nSigGenes < nRefCells) {
    stop(
      "There are only ", nSigGenes, " signature genes",
      " matching common genes between bulk and reference profiles,",
      " but there should be more signature genes than reference cells"
    )
  }

  if (scaleExprs) {
    if (length(commonGenes) < 2e3) {
      warning(
        "there are few genes in common between the bulk samples and ",
        "reference cells:", length(commonGenes), ", so the data scaling ",
        "might be an issue"
      )
    }
    # The value of 2e3 is arbitrary, but should be a respectable number for the
    # data renormalization.
    bulk <- scaleCounts(bulk, sigGenes, commonGenes)$counts
    temp <- scaleCounts(refProfiles, sigGenes, commonGenes)
    refProfiles <- temp$counts
    if (with_w) {
      refProfiles.var <- scaleCounts(refProfiles.var, sigGenes,
        normFact = temp$normFact
      )$counts
    }
    # the refProfiles.var is normalized by the same factors as refProfiles.
  } else {
    bulk <- bulk[sigGenes, ]
    refProfiles <- refProfiles[sigGenes, ]
    if (with_w) {
      refProfiles.var <- refProfiles.var[sigGenes, ]
    }
  }

  if (is.null(mRNA_cell)) {
    mRNA_cell <- IOBR::mRNA_cell_default
  }

  if (!is.null(mRNA_cell_sub)) {
    if (is.null(names(mRNA_cell_sub)) || !is.numeric(mRNA_cell_sub)) {
      stop("When mRNA_cell_sub is given, it needs to be a named numeric vector")
    }
    mRNA_cell[names(mRNA_cell_sub)] <- mRNA_cell_sub
  }

  minFun <- function(x, A, b, w) {
    # Basic minimization function used to minimize the squared sum of the error
    # between the fit and observed value (A*x - b). We also give a weight, w,
    # for each gene to give more or less importance to the fit of each (can use
    # a value 1 if don't want to give any weight).
    return(sum((w * (A %*% x - b)^2), na.rm = TRUE))
  }
  minFun.range <- function(x, A, b, A.var) {
    # Other minimization function where we don't use weights but instead give
    # also the variability on A as input and we'll compute for each gene the
    # min / max value of the pred based on this variability to keep the smallest
    # value. I.e. we compute a range of the values that y_pred can take for each
    # gene based on the ref profile variability and when the y_true is within
    # this range then we return an error of 0 for this gene and if y_true is
    # outside of the range, we keep the "smallest error" (i.e. value most
    # nearby of the possible range).
    val.max <- (A + A.var) %*% x - b
    val.min <- (A - A.var) %*% x - b
    cErr <- rep(0, length(b))
    outOfRange <- (sign(val.max) * sign(val.min) == 1)
    cErr[outOfRange] <- pmin(abs(val.max[outOfRange]), abs(val.min[outOfRange]))
    return(sum(cErr, na.rm = TRUE))
  }

  if (with_w && !rangeBasedOptim) {
    # Computing the weight to give to each gene
    w <- rowSums(refProfiles / (refProfiles.var + 1e-12), na.rm = TRUE)
    # 1e-12 to avoid divisions by 0: like this if refProfiles and refProfiles.var
    # both are 0 for a given element, it will result in a weight of 0.
    med_w <- stats::median(w[w > 0], na.rm = TRUE)
    w[w > 100 * med_w] <- 100 * med_w
    # Set a limit for the big w to still not give too much weights to these genes
  } else {
    w <- 1
  }

  # Defining the constraints for the fit of the proportions.
  if (withOtherCells) {
    cMin <- 0
  } else {
    cMin <- 0.99
    # So that when we assume no cells without known reference profile are
    # present, we ask the sum of the mRNA proportions of known cells to be at
    # least 0.99.
  }
  cMax <- 1
  ui <- diag(nRefCells)
  ci <- rep(0, nRefCells)
  if (constrainedSum) {
    ui <- rbind(ui, rep(1, nRefCells), rep(-1, nRefCells))
    ci <- c(ci, cMin, -cMax)
  }
  # ui and ci define the constraints, in the form "ui %*% x - ci >= 0".
  # The first constraints are that the proportion of each cell must be
  # positive; then if constrainedSum is true, we want the sum of all proportions
  # to be bigger than cMin and this sum must also be smaller or equal to cMax.
  cInitProp <- (min(1, cMax) - 1e-5) / nRefCells
  # used as an initial guess for the optimized - start with all cell fractions
  # equal. We added the -1e-5 in cInitProp because the optimizer needs to have
  # the initial guess inside the admissible region and not on its boundary

  # Estimating for each sample the proportion of the mRNA per cell type.
  tempPropPred <- lapply(1:nSamples, FUN = function(cSample) {
    b <- bulk[, cSample]
    if (!rangeBasedOptim) {
      fit <- stats::constrOptim(
        theta = rep(cInitProp, nRefCells), f = minFun,
        grad = NULL, ui = ui, ci = ci, A = refProfiles, b = b, w = w
      )
    } else {
      fit <- stats::constrOptim(
        theta = rep(cInitProp, nRefCells), f = minFun.range,
        grad = NULL, ui = ui, ci = ci, A = refProfiles, b = b, A.var = refProfiles.var
      )
    }
    fit$x <- fit$par
    if (!withOtherCells) {
      fit$x <- fit$x / sum(fit$x, na.rm = TRUE)
    }
    # So that the sum is really equal to 1 even if the best pred was giving
    # slightly lower values when we force the system to have only known cells.

    # Checking how well the estimated proportions predict the gene expression
    b_estimated <- refProfiles %*% fit$x

    if (nSigGenes > 2) {
      suppressWarnings(corSp.test <- stats::cor.test(b, b_estimated, method = "spearman"))
      # Use suppressWarnings to avoid the warnings that p-value isn't exact when
      # ties are present in the data.
      corPear.test <- stats::cor.test(b, b_estimated, method = "pearson")
    } else {
      # cannot compute correlations with less than 3 observations.
      corSp.test <- corPear.test <- list()
      corSp.test$estimate <- corSp.test$p.value <- corPear.test$estimate <-
        corPear.test$p.value <- NA
    }
    regLine <- stats::lm(b_estimated ~ b)
    regLine_through0 <- stats::lm(b_estimated ~ b + 0)
    if (!rangeBasedOptim) {
      rmse_pred <- sqrt(minFun(x = fit$x, A = refProfiles, b = b, w = w) / nSigGenes)
      rmse_0 <- sqrt(minFun(x = rep(0, nRefCells), A = refProfiles, b = b, w = w) / nSigGenes)
    } else {
      rmse_pred <- sqrt(minFun.range(
        x = fit$x, A = refProfiles, b = b,
        A.var = refProfiles.var
      ) / nSigGenes)
      rmse_0 <- sqrt(minFun.range(
        x = rep(0, nRefCells), A = refProfiles, b = b,
        A.var = refProfiles.var
      ) / nSigGenes)
    }
    gof <- data.frame(fit$convergence, ifelse(is.null(fit$message), "", fit$message),
      rmse_pred, rmse_0,
      # to only have sum((w*b)^2) or corresponding value based
      # on the minFun, in the worst case possible.
      corSp.test$estimate, corSp.test$p.value,
      corPear.test$estimate, corPear.test$p.value,
      regLine$coefficients[2], regLine$coefficients[1],
      regLine_through0$coefficients[1], sum(fit$x),
      stringsAsFactors = FALSE
    )

    return(list(mRNAProportions = fit$x, fit.gof = gof))
  })
  mRNAProportions <- do.call(rbind, lapply(
    tempPropPred,
    function(x) x$mRNAProportions
  ))
  dimnames(mRNAProportions) <- list(samplesNames, refCellsNames)

  fit.gof <- do.call(rbind, lapply(tempPropPred, function(x) x$fit.gof))
  dimnames(fit.gof) <- list(samplesNames, c(
    "convergeCode", "convergeMessage", "RMSE_weighted",
    "Root_mean_squared_geneExpr_weighted",
    "spearmanR", "spearmanP", "pearsonR", "pearsonP", "regline_a_x",
    "regline_b", "regline_a_x_through0", "sum_mRNAProportions"
  ))
  # Some matrix giving information on the goodness of the fit.

  if (any(fit.gof$convergeCode != 0)) {
    warning(
      "The optimization didn't fully converge for some samples:\n",
      paste(samplesNames[fit.gof$convergeCode != 0], collapse = "; "),
      "\n - check fit.gof for the convergeCode and convergeMessage"
    )
  }

  if (withOtherCells) {
    mRNAProportions <- cbind(mRNAProportions, otherCells = 1 - rowSums(mRNAProportions))
  }
  # Adding a row to the proportion matrix corresponding to the other /
  # uncharacterized / cancer cells.

  tInds <- match(colnames(mRNAProportions), names(mRNA_cell))
  if (anyNA(tInds)) {
    defaultInd <- match("default", names(mRNA_cell))
    if (is.na(defaultInd)) {
      tStr <- paste(
        " and no default value is given for this mRNA per cell,",
        "so we cannot estimate the cellFractions, only",
        "the mRNA proportions"
      )
    } else {
      tStr <- paste(
        " - using the default value of", mRNA_cell[defaultInd],
        "for these but this might bias the true cell proportions from",
        "all cell types."
      )
    }
    warning(
      "mRNA_cell value unknown for some cell types: ",
      paste(colnames(mRNAProportions)[is.na(tInds)], collapse = ", "),
      tStr
    )
    tInds[is.na(tInds)] <- defaultInd
  }
  cellFractions <- t(t(mRNAProportions) / mRNA_cell[tInds])
  cellFractions <- cellFractions / rowSums(cellFractions, na.rm = FALSE)
  # So that the cell fractions sum up to 1. In case some values were NA (due
  # to unknown mRNA_cell value for some cell types without default value), then
  # the full cellFractions will be NA - if we used na.rm=T, then we'd have the
  # sum of the other than "NA" cells equal to 1 as if this "NA" cell was not
  # present in the sample.

  return(list(
    mRNAProportions = mRNAProportions, cellFractions = cellFractions,
    fit.gof = fit.gof
  ))
}

#' Scaling raw counts from each sample.
#'
#' Normalizing the sum of counts from each sample to 1e6.
#'
#' Function taking a matrix (\emph{genes} x \emph{samples}) of counts as input
#' and returning the scaled counts for the subset signature genes (or all genes
#' if it sigGenes is \code{NULL}), with the scaled counts computed based on all
#' the 'renormGenes' (or all genes if it is NULL). The renormalization is made
#' independently for each sample, to have the sum of each columns over the
#' renormGenes equal to 1e6.
#'
#' normFact, if not null, is used as the normalization factor instead of
#'  the colSums (used to renormalize the refProfiles.var by the same amount
#'  than the refProfiles).
#'
#' @keywords internal
scaleCounts <- function(counts, sigGenes = NULL, renormGenes = NULL, normFact = NULL) {
  if (is.null(sigGenes)) {
    sigGenes <- 1:nrow(counts)
  }

  if (is.null(normFact)) {
    if (is.null(renormGenes)) {
      renormGenes <- 1:nrow(counts)
    }
    normFact <- colSums(counts[renormGenes, , drop = FALSE], na.rm = TRUE)
  }
  counts <- t(t(counts[sigGenes, , drop = FALSE]) / normFact) * 1e6
  # Need to take the transpose so that the division is made on the correct
  # elements
  return(list(counts = counts, normFact = normFact))
}

#' Merging the duplicates from the input matrix.
#'
#' In case there are some duplicate rownames in the input matrix (\emph{mat}),
#' this function will return a similar matrix with unique rownames and the
#' duplicate cases will be merged together based on the median values. When
#' \emph{warn} is true a warning message will be written in case there are
#' duplicates, this warning message will include the \emph{in_type} string to
#' indicate the current type of the input matrix.
#'
#' @param mat matrix
#' @param warn warning
#' @param in_type default is null
#'
#' @keywords internal
merge_duplicates <- function(mat, warn = TRUE, in_type = NULL) {
  dupl <- duplicated(rownames(mat))
  if (sum(dupl) > 0) {
    dupl_genes <- unique(rownames(mat)[dupl])
    if (warn) {
      warning(
        "There are ", length(dupl_genes), " duplicated gene names",
        ifelse(!is.null(in_type), paste(" in the", in_type), ""),
        ". We'll use the median value for each of these cases."
      )
    }
    mat_dupl <- mat[rownames(mat) %in% dupl_genes, , drop = F]
    mat_dupl_names <- rownames(mat_dupl)
    mat <- mat[!dupl, , drop = F]
    # First put the dupl cases in a separate matrix and keep only the unique
    # gene names in the mat matrix.
    mat[dupl_genes, ] <- t(sapply(dupl_genes, FUN = function(cgene) {
      apply(mat_dupl[mat_dupl_names == cgene, , drop = F], MARGIN = 2, FUN = median)
    }))
  }
  return(mat)
}
