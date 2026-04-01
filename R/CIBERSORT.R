# CIBERSORT R script v1.04 (last updated 10-24-2016)
# Note: Signature matrix construction is not currently available;
#   use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#   R v3.0 or later.
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt

#' Core Algorithm for CIBERSORT Deconvolution
#'
#' @description
#' Performs nu-regression using support vector machines (SVM) to estimate
#' cell type proportions. This is the core computational engine of CIBERSORT,
#' using nu-SVR with linear kernel to decompose mixed gene expression signals.
#'
#' @param X Matrix or data frame containing signature matrix (predictor variables).
#'   Rows are genes, columns are cell types.
#' @param y Numeric vector containing the mixture sample expression (response variable).
#' @param absolute Logical indicating whether to use absolute space for weights.
#'   Default is FALSE (relative proportions).
#' @param abs_method String specifying the method for absolute space weights:
#'   `"sig.score"` or `"no.sumto1"`.
#'
#' @return List containing:
#' \describe{
#'   \item{w}{Estimated cell type weights/proportions}
#'   \item{mix_rmse}{Root mean squared error of the deconvolution}
#'   \item{mix_r}{Correlation coefficient between observed and predicted mixture}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' X <- matrix(rnorm(100), nrow = 10)
#' y <- rnorm(10)
#' result <- CoreAlg(X, y, absolute = FALSE, abs_method = "sig.score")
#' }
CoreAlg <- function(X, y, absolute, abs_method) {
  rlang::check_installed("e1071")

  svn_itor <- 3
  nus <- c(0.25, 0.5, 0.75)

  res <- function(i) {
    model <- e1071::svm(X, y,
      type = "nu-regression", kernel = "linear",
      nu = nus[i], scale = FALSE
    )
    model
  }

  if (Sys.info()["sysname"] == "Windows") {
    out <- parallel::mclapply(1:svn_itor, res, mc.cores = 1)
  } else {
    out <- parallel::mclapply(1:svn_itor, res, mc.cores = 1)
  }

  nusvm <- numeric(svn_itor)
  corrv <- numeric(svn_itor)

  for (t in seq_len(svn_itor)) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt(mean((k - y)^2))
    corrv[t] <- cor(k, y)
  }

  mn <- which.min(nusvm)
  model <- out[[mn]]

  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0

  if (!absolute || abs_method == "sig.score") {
    w <- q / sum(q)
  } else if (absolute && abs_method == "no.sumto1") {
    w <- q
  }

  list(w = w, mix_rmse = nusvm[mn], mix_r = corrv[mn])
}

#' Permutation Test for CIBERSORT
#'
#' @description
#' Performs permutation-based sampling to generate an empirical null distribution
#' of correlation coefficients for p-value calculation in CIBERSORT analysis.
#' Randomly samples from the mixture data to create null distributions.
#'
#' @param perm Integer. Number of permutations to perform (≥100 recommended for
#'   reliable p-value estimation).
#' @param X Matrix or data frame containing signature matrix (predictor variables).
#' @param Y Numeric vector containing the mixture sample expression.
#' @param absolute Logical indicating whether to use absolute space for weights.
#' @param abs_method String specifying the method for absolute space weights:
#'   `"sig.score"` or `"no.sumto1"`.
#' @param seed Integer. Random seed for reproducibility. If NULL (default),
#'   uses current random state.
#'
#' @return List containing:
#' \describe{
#'   \item{dist}{Numeric vector of correlation coefficients from permutations}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' X <- matrix(rnorm(100), nrow = 10)
#' Y <- rnorm(10)
#' result <- doPerm(1000, X, Y, absolute = FALSE, abs_method = "sig.score")
#' }
doPerm <- function(perm, X, Y, absolute, abs_method, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (perm < 1) {
    cli::cli_abort("{.arg perm} must be at least 1, got {perm}")
  }

  Ylist <- as.list(data.matrix(Y))
  dist <- numeric(perm)

  for (itor in seq_len(perm)) {
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])
    yr <- (yr - mean(yr)) / sd(yr)
    result <- CoreAlg(X, yr, absolute, abs_method)
    dist[itor] <- result$mix_r
  }

  list(dist = dist)
}


#' Parallel Permutation Test for CIBERSORT
#'
#' @description
#' Parallel version of doPerm. Performs permutation-based sampling and runs
#' the CoreAlg function iteratively using multiple CPU cores to accelerate
#' computation. This function generates an empirical null distribution of
#' correlation coefficients for p-value calculation in CIBERSORT analysis.
#'
#' @param perm1 Integer. Number of permutations to perform (≥100 recommended).
#' @param X1 Matrix or data frame. Signature matrix (cell type GEP barcode).
#' @param Y1 Matrix. Mixture file containing gene expression profiles.
#' @param absolute1 Logical. Whether to run in absolute mode (default: FALSE).
#' @param abs_method1 String. Method for absolute mode: 'sig.score' or 'no.sumto1'.
#' @param num_cores1 Integer. Number of CPU cores for parallel computation (default: 2).
#' @param seed Integer. Random seed for reproducibility. If NULL (default),
#'   uses current random state. Set to a specific value (e.g., 123) for
#'   reproducible results across runs.
#'
#' @return List containing:
#' \describe{
#'   \item{dist}{Numeric vector of correlation coefficients from permutations,
#'     representing the empirical null distribution.}
#' }
#'
#' @details
#' This function utilizes the \code{foreach} and \code{doParallel} packages
#' to distribute permutation iterations across multiple cores. It automatically
#' handles cluster setup/teardown via \code{on.exit()} to prevent resource leaks.
#'
#' Note: Windows users may experience slower performance due to socket-based
#' parallelization (PSOCK) versus forking on Unix systems.
#'
#' @seealso \code{\link{doPerm}} for the sequential version, \code{\link{CoreAlg}}, \code{\link{CIBERSORT}}
#'
#' @importFrom parallel makeCluster stopCluster
#' @export
#'
#' @examples
#' \donttest{
#' X <- matrix(rnorm(1000), nrow = 100)
#' Y <- matrix(rnorm(500), nrow = 100)
#' rownames(X) <- rownames(Y) <- paste0("Gene", 1:100)
#'
#' result <- parallel_doperm(
#'   perm1 = 100, X1 = X, Y1 = Y,
#'   absolute1 = FALSE, abs_method1 = "sig.score", num_cores1 = 2
#' )
#' str(result$dist)
#' }
parallel_doperm <- function(perm1, X1, Y1, absolute1, abs_method1,
                            num_cores1 = 2, seed = NULL) {
  rlang::check_installed("foreach", reason = "for parallel permutation")
  rlang::check_installed("doParallel", reason = "for parallel backend")

  if (perm1 < 1) {
    cli::cli_abort("{.arg perm1} must be at least 1, got {perm1}")
  }
  if (num_cores1 < 1) {
    cli::cli_abort("{.arg num_cores1} must be at least 1, got {num_cores1}")
  }

  `%dopar%` <- foreach::`%dopar%`

  if (!is.null(seed)) set.seed(seed)
  n <- dim(X1)[1]
  Ylist_len <- length(as.list(data.matrix(Y1)))
  perm_indices <- lapply(seq_len(perm1), function(i) sample(Ylist_len, n))

  cl <- parallel::makeCluster(num_cores1)
  on.exit(
    {
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    },
    add = TRUE
  )

  doParallel::registerDoParallel(cl)

  dist <- foreach::foreach(
    itor = seq_len(perm1),
    .combine = rbind,
    .export = c("perm_indices", "Y1", "X1", "absolute1", "abs_method1"),
    .packages = "e1071"
  ) %dopar% {
    Ylist <- as.list(data.matrix(Y1))
    yr <- as.numeric(Ylist[perm_indices[[itor]]])
    yr <- (yr - mean(yr)) / sd(yr)

    result <- CoreAlg(X1, yr, absolute1, abs_method1)
    result$mix_r
  }

  list(dist = dist)
}


#' CIBERSORT Deconvolution Algorithm
#'
#' An analytical tool to estimate cell type abundances in mixed cell populations
#' using gene expression data.

#' @param sig_matrix Cell type GEP barcode matrix: row 1 = sample labels;
#'   column 1 = gene symbols; no missing values; default = LM22.txt download
#'   from CIBERSORT (https://cibersort.stanford.edu/runcibersort.php)
#' @param mixture_file GEP matrix: row 1 = sample labels; column 1 = gene
#'   symbols; no missing values
#' @param perm Set permutations for statistical analysis
#'   (>=100 permutations recommended).
#' @param QN Quantile normalization of input mixture (default = TRUE)
#' @param absolute Run CIBERSORT in absolute mode (default = FALSE)
#'   - note that cell subsets will be scaled by their absolute levels and
#'     will not be represented as fractions (to derive the default output,
#'     normalize absolute levels such that they sum to 1 for each mixture sample)
#'   - the sum of all cell subsets in each mixture sample will be added to the
#'     output ('Absolute score'). If LM22 is used, this score will capture
#'     total immune content.
#' @param abs_method If absolute is set to TRUE, choose method:
#'   'no.sumto1' or 'sig.score'
#'   - sig.score = for each mixture sample, define S as the median expression
#'     level of all genes in the signature matrix divided by the median
#'     expression level of all genes in the mixture. Multiple cell subset
#'     fractions by S.
#'   - no.sumto1 = remove sum to 1 constraint
#' @param parallel Logical. Enable parallel execution? (default = FALSE)
#' @param num_cores Integer. Number of cores to use when \code{parallel = TRUE}
#'   (default = 2)
#' @param seed Integer. Random seed for reproducible permutation testing.
#'   If \code{NULL} (default), uses current random state. Set to a specific
#'   value (e.g., 123) for reproducible results across runs. Applies to both
#'   parallel and serial permutation.
#' @author Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' @return A matrix object containing the estimated cibersort-cell fractions,
#'   p-values, correlation coefficients, and RMSE values.
#' @export
#' @import parallel
#' @import preprocessCore
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import purrr
#' @import stringr
#' @examples
#' \donttest{
#' # Create simulated data matching LM22 signature matrix gene names
#' data(lm22)
#' common_genes <- rownames(lm22)[1:500]
#' sim_mixture <- as.data.frame(matrix(
#'   rnorm(length(common_genes) * 10, mean = 5, sd = 2),
#'   nrow = length(common_genes), ncol = 10
#' ))
#' rownames(sim_mixture) <- common_genes
#' colnames(sim_mixture) <- paste0("Sample", 1:10)
#' result <- CIBERSORT(
#'   sig_matrix = lm22,
#'   mixture_file = sim_mixture,
#'   perm = 10, QN = FALSE, absolute = FALSE,
#'   parallel = FALSE
#' )
#' head(result)
#' }
CIBERSORT <- function(sig_matrix = NULL, mixture_file, perm, QN = TRUE,
                      absolute = FALSE, abs_method = "sig.score",
                      parallel = FALSE, num_cores = 2, seed = NULL) {
  if (is.null(sig_matrix)) {
    sig_matrix <- load_data("lm22")
  }
  # Input validation
  if (length(intersect(rownames(mixture_file), rownames(sig_matrix))) == 0) {
    stop("None identical gene between eset and reference had been found.
         Check your eset using: intersect(rownames(eset), rownames(reference))")
  }

  if (absolute && abs_method != "no.sumto1" && abs_method != "sig.score") {
    stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  }

  # Read in data
  X <- sig_matrix
  Y <- rownames_to_column(mixture_file, var = "symbol")

  # Handle duplicated gene symbols
  dups <- dim(Y)[1] - length(unique(Y[, 1]))
  if (dups > 0) {
    warning(paste(dups, "duplicated gene symbol(s) found in mixture file!"))
    rownames(Y) <- make.names(Y[, 1], unique = TRUE)
  } else {
    rownames(Y) <- Y[, 1]
  }
  Y <- Y[, -1, drop = FALSE]

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  # Order
  X <- X[order(rownames(X)), , drop = FALSE]
  Y <- Y[order(rownames(Y)), , drop = FALSE]

  P <- perm

  # Anti-log if max < 50 in mixture file
  if (max(Y) < 50) {
    Y <- 2^Y
  }

  # Quantile normalization of mixture file
  if (QN == TRUE) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Yq <- normalize.quantiles(Y)
    Y <- matrix(Yq,
      nrow = length(tmpr), ncol = length(tmpc),
      dimnames = list(tmpr, tmpc)
    )
  }

  # Store original mixtures
  Yorig <- Y
  Ymedian <- max(median(Yorig), 1)

  # Intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, , drop = FALSE]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY, , drop = FALSE]

  # Sanity check
  if (length(dim(Y)) != 2) {
    stop(
      "`mixture_file` became non-2D after preprocessing. ",
      "This usually occurs when it has only one column and subsetting ",
      "dropped dimensions. ",
      "Please update package or ensure `mixture_file` is a matrix with drop = FALSE."
    )
  }

  # Standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  # Parallel setup with safety checks
  use_parallel_mixtures <- parallel

  if (parallel) {
    rlang::check_installed("BiocParallel", reason = "for parallel mode")

    # Cap cores to prevent system overload
    max_available <- parallel::detectCores()
    if (is.na(max_available)) max_available <- 2
    safe_cores <- max(1, max_available - 1)

    if (num_cores > safe_cores) {
      warning(
        "Requested ", num_cores, " cores but only ", safe_cores,
        " available/safe. Reducing to ", safe_cores, "."
      )
      num_cores <- safe_cores
    }

    BPPARAM <- BiocParallel::MulticoreParam(workers = num_cores)
  }

  # Empirical null distribution of correlation coefficients
  if (P > 0) {
    if (parallel && exists("parallel_doperm", mode = "function")) {
      # Use parallel permutation (consumes parallel resources)
      kk <- parallel_doperm(
        perm1 = P, X1 = X, Y1 = Y,
        absolute1 = absolute, abs_method1 = abs_method,
        num_cores1 = num_cores, seed = seed
      )
      nulldist <- sort(kk$dist)
      # CRITICAL: Disable parallel for mixtures to avoid nested parallelism
      use_parallel_mixtures <- FALSE
    } else {
      # Serial permutation
      nulldist <- sort(doPerm(P, X, Y, absolute, abs_method, seed = seed)$dist)
    }
  } else {
    nulldist <- NULL
  }

  header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
  if (absolute) {
    header <- c(header, paste0("Absolute score (", abs_method, ")"))
  }

  # Core job function with explicit parameters (Windows-safe)
  coreJob <- function(itor, Y, X, absolute, abs_method, Ymedian, P, nulldist) {
    y <- Y[, itor, drop = TRUE]
    y <- (y - mean(y)) / sd(y)
    result <- CoreAlg(X, y, absolute, abs_method)
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    if (absolute && abs_method == "sig.score") {
      w <- w * median(Y[, itor]) / Ymedian
    }
    pval <- if (P > 0) {
      1 - which.min(abs(nulldist - mix_r)) / length(nulldist)
    } else {
      9999
    }

    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    if (absolute) out <- c(out, sum(w))
    out
  }

  mixtures <- ncol(Y)

  # Process mixtures
  if (use_parallel_mixtures && mixtures > 1) {
    # Parallel processing (Windows-safe via explicit variable passing)
    output_list <- BiocParallel::bplapply(
      seq_len(mixtures),
      function(i) {
        coreJob(i,
          Y = Y, X = X, absolute = absolute,
          abs_method = abs_method, Ymedian = Ymedian,
          P = P, nulldist = nulldist
        )
      },
      BPPARAM = BPPARAM
    )
    output <- do.call(rbind, output_list)
  } else {
    # Serial processing (safe fallback)
    output_list <- lapply(
      seq_len(mixtures),
      function(i) {
        coreJob(i,
          Y = Y, X = X, absolute = absolute,
          abs_method = abs_method, Ymedian = Ymedian,
          P = P, nulldist = nulldist
        )
      }
    )
    output <- do.call(rbind, output_list)
  }

  # Return matrix object containing all results
  obj <- rbind(header, output)
  obj <- obj[, -1, drop = FALSE]
  obj <- obj[-1, , drop = FALSE]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  if (!absolute) {
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  } else {
    colnames(obj) <- c(
      colnames(X), "P-value", "Correlation", "RMSE",
      paste0("Absolute score (", abs_method, ")")
    )
  }
  obj
}
