# CIBERSORT R script v1.04 (last updated 10-24-2016)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.packges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)
#
#       Options:
#       i)   perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii)  QN = Quantile normalization of input mixture (default = TRUE)
#       iii) absolute = Run CIBERSORT in absolute mode (default = FALSE)
#               - note that cell subsets will be scaled by their absolute levels and will not be
#                 represented as fractions (to derive the default output, normalize absolute
#                 levels such that they sum to 1 for each mixture sample)
#               - the sum of all cell subsets in each mixture sample will be added to the ouput
#                 ('Absolute score'). If LM22 is used, this score will capture total immune content.
#       iv)  abs_method = if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#               - sig.score = for each mixture sample, define S as the median expression
#                 level of all genes in the signature matrix divided by the median expression
#                 level of all genes in the mixture. Multiple cell subset fractions by S.
#               - no.sumto1 = remove sum to 1 constraint
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


# dependencies
# library(e1071)
# library(parallel)
# library(preprocessCore)

# Core algorithm
#' Perform Nu-Regression Using Support Vector Machines
#'
#' This function performs nu-regression using support vector machines (SVM) and calculates weights,
#' root mean squared error (RMSE), and correlation coefficient (R).
#'
#' @param X A matrix or data frame containing the predictor variables.
#' @param y A numeric vector containing the response variable.
#' @param absolute Logical indicating whether to use absolute space for weights. Default is FALSE.
#' @param abs_method String specifying the method for absolute space weights: "sig.score" or "no.sumto1".
#'
#' @return A list containing the weights (`w`), root mean squared error (`mix_rmse`), and correlation coefficient (`mix_r`).
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10)
#' y <- rnorm(10)
#' result <- CoreAlg(X, y, absolute = FALSE, abs_method = "sig.score")
#'
CoreAlg <- function(X, y, absolute, abs_method) {
  rlang::check_installed("e1071")
  # try different values of nu
  svn_itor <- 3

  res <- function(i) {
    if (i == 1) {
      nus <- 0.25
    }
    if (i == 2) {
      nus <- 0.5
    }
    if (i == 3) {
      nus <- 0.75
    }
    model <- e1071::svm(X, y, type = "nu-regression", kernel = "linear", nu = nus, scale = F)
    model
  }

  if (Sys.info()["sysname"] == "Windows") {
    out <- mclapply(1:svn_itor, res, mc.cores = 1)
  } else {
    out <- mclapply(1:svn_itor, res, mc.cores = svn_itor)
  }

  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)

  # do cibersort
  t <- 1
  while (t <= svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  # pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  # get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0
  if (!absolute || abs_method == "sig.score") w <- (q / sum(q)) # relative space (returns fractions)
  if (absolute && abs_method == "no.sumto1") w <- q # absolute space (returns scores)

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}

# do permutations
#' Title
#' @description doPerm performs permutation-based sampling and runs the CoreAlg function iteratively.
#' @param perm Number of permutations to perform.
#' @param X Input matrix or data frame containing the predictor variables.
#' @param Y Numeric vector containing the response variable.
#' @param absolute Logical value indicating whether to use absolute space or relative space for the weights.
#' @param abs_method String indicating the method to calculate the weights in absolute space. Can be either 'sig.score' or 'no.sumto1'.
#' @param seed Integer. Random seed for reproducibility. If NULL (default), 
#'   uses current random state. Set to a specific value for reproducible results.
#'   
#' @return A list containing the distribution of correlation coefficients from the permutations.
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10)
#' y <- rnorm(10)
#' result <- doPerm(1000, X, Y, absolute = FALSE, abs_method = "sig.score")
#'
doPerm <- function(perm, X, Y, absolute, abs_method,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while (itor <= perm) {
    # print(itor)

    # random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])

    # standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    # run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)

    mix_r <- result$mix_r

    # store correlation
    if (itor == 1) {
      dist <- mix_r
    } else {
      dist <- rbind(dist, mix_r)
    }

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
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
#' @return 
#' A list containing:
#' \item{dist}{Numeric vector of correlation coefficients from permutations, 
#'   representing the empirical null distribution.}
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
#' \dontrun{
#' # Prepare data (using example data structure from CIBERSORT)
#' X <- matrix(rnorm(1000), nrow = 100)  # Signature matrix
#' Y <- matrix(rnorm(500), nrow = 100)   # Mixture matrix
#' rownames(X) <- rownames(Y) <- paste0("Gene", 1:100)
#' 
#' # Run parallel permutation (using 2 cores)
#' result <- parallel_doperm(
#'   perm1 = 100, 
#'   X1 = X, 
#'   Y1 = Y, 
#'   absolute1 = FALSE, 
#'   abs_method1 = "sig.score",
#'   num_cores1 = 2
#' )
#' str(result$dist)
#' }
parallel_doperm <- function(perm1, X1, Y1, absolute1, abs_method1, num_cores1 = 2,
                            seed=NULL) {
  # 检查依赖包
  rlang::check_installed("foreach", reason = "for parallel permutation")
  rlang::check_installed("doParallel", reason = "for parallel backend")
  
  `%dopar%` <- foreach::`%dopar%`
  # 预生成随机索引（确定性）
  if (!is.null(seed)) set.seed(seed)
  n <- dim(X1)[1]
  Ylist_len <- length(as.list(data.matrix(Y1)))
  perm_indices <- lapply(1:perm1, function(i) sample(Ylist_len, n))
  
 
  # 创建集群并使用 on.exit 确保资源释放
  # 注意：makeCluster 是 parallel 包的函数，不是 doParallel
  cl <- parallel::makeCluster(num_cores1)
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()  # 重置为顺序模式，避免影响后续代码
  }, add = TRUE)
  
  doParallel::registerDoParallel(cl)
  
  # 并行执行置换检验（使用预计算索引，无随机性）
  dist <- foreach::foreach(
    itor = 1:perm1,
    .combine = rbind,
    .export = c("perm_indices", "Y1", "X1", "absolute1", "abs_method1"),  # 必须导出 perm_indices
    .packages = 'e1071'
  ) %dopar% {
    # 使用预计算的确定性索引
    Ylist <- as.list(data.matrix(Y1))
    yr <- as.numeric(Ylist[perm_indices[[itor]]])  
    yr <- (yr - mean(yr)) / sd(yr)
    
    result <- CoreAlg(X1, yr, absolute1, abs_method1)
    result$mix_r
  }
  
  list(dist = dist)
}



#' CIBERSORT is an analytical tool developed by Newman et al. to provide an estimation of the abundances of member cell types in a mixed cell population, using gene expression data.

#' @param sig_matrix  Cell type GEP barcode matrix: row 1 = sample labels; column 1 = gene symbols; no missing values; default =LM22.txt download from CIBERSORT (https://cibersort.stanford.edu/runcibersort.php)
#' @param mixture_file  GEP matrix: row 1 = sample labels; column 1 = gene symbols; no missing values
#' @param perm Set permutations for statistical analysis (≥100 permutations recommended).
#' @param QN Quantile normalization of input mixture (default = TRUE)
#' @param absolute  Run CIBERSORT in absolute mode (default = FALSE)
#' - note that cell subsets will be scaled by their absolute levels and will not be represented as fractions
#' (to derive the default output, normalize absolute levels such that they sum to 1 for each mixture sample)
#' - the sum of all cell subsets in each mixture sample will be added to the ouput ('Absolute score').
#' If LM22 is used, this score will capture total immune content.
#' @param abs_method  if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#' - sig.score = for each mixture sample, define S as the median expression
#' level of all genes in the signature matrix divided by the median expression
#' level of all genes in the mixture. Multiple cell subset fractions by S.
#' - no.sumto1 = remove sum to 1 constraint
#' @param parallel Logical. Enable parallel execution? (default = FALSE)
#' @param num_cores Integer. Number of cores to use when \code{parallel = TRUE} (default = 2)
#' @param seed Integer. Random seed for reproducible permutation testing. 
#'   If \code{NULL} (default), uses current random state. Set to a specific 
#'   value (e.g., 123) for reproducible results across runs. Applies to both 
#'   parallel and serial permutation.
#' @author Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' @return A matrix object containing the estimated cibersort-cell fractions, p-values, correlation coefficients, and RMSE values.
#' @export
#' @import parallel
#' @import preprocessCore
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import purrr
#' @import stringr
#' @examples
#' data("eset_gse62254", package = "IOBR")
#' cibersort <- CIBERSORT(sig_matrix = lm22,
#'   mixture_file = eset_gse62254, perm = 100, 
#'   QN = TRUE, absolute = FALSE)
#' head(cibersort)
CIBERSORT <- function(sig_matrix = lm22, mixture_file, perm, QN = TRUE, absolute,
                      abs_method = "sig.score", parallel = FALSE, num_cores = 2,
                      seed=NULL){
  
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
    warning(paste(dups, " duplicated gene symbol(s) found in mixture file!", sep = ""))
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
    Y <- matrix(Yq, nrow = length(tmpr), ncol = length(tmpc),
                dimnames = list(tmpr, tmpc))
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
    stop("`mixture_file` became non-2D after preprocessing. ",
         "This usually occurs when it has only one column and subsetting dropped dimensions. ",
         "Please update package or ensure `mixture_file` is a matrix with drop = FALSE.")
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
      warning("Requested ", num_cores, " cores but only ", safe_cores, 
              " available/safe. Reducing to ", safe_cores, ".")
      num_cores <- safe_cores
    }
    
    BPPARAM <- BiocParallel::MulticoreParam(workers = num_cores)
  }
  
  # Empirical null distribution of correlation coefficients
  if (P > 0) {
    if (parallel && exists("parallel_doperm", mode = "function")) {
      # Use parallel permutation (consumes parallel resources)
      kk <- parallel_doperm(perm1 = P, X1 = X, Y1 = Y,
                            absolute1 = absolute, abs_method1 = abs_method,
                            num_cores1 = num_cores,seed=seed)   
      nulldist <- sort(kk$dist)
      # CRITICAL: Disable parallel for mixtures to avoid nested parallelism
      use_parallel_mixtures <- FALSE
    } else {
      # Serial permutation
      nulldist <- sort(doPerm(P, X, Y, absolute, abs_method,seed = seed)$dist)
    }
  } else {
    nulldist <- NULL
  }
  
  header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
  if (absolute) header <- c(header, paste("Absolute score (", abs_method, ")", sep = ""))
  
  # Core job function with explicit parameters (Windows-safe)
  coreJob <- function(itor, Y, X, absolute, abs_method, Ymedian, P, nulldist) {
    y <- Y[, itor, drop = TRUE]
    y <- (y - mean(y)) / sd(y)
    result <- CoreAlg(X, y, absolute, abs_method)
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    if (absolute && abs_method == "sig.score")
      w <- w * median(Y[, itor]) / Ymedian
    pval <- if (P > 0)
      1 - which.min(abs(nulldist - mix_r)) / length(nulldist) else 9999
    
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
      function(i) coreJob(i, Y = Y, X = X, absolute = absolute, 
                          abs_method = abs_method, Ymedian = Ymedian, 
                          P = P, nulldist = nulldist),
      BPPARAM = BPPARAM
    )
    output <- do.call(rbind, output_list)
  } else {
    # Serial processing (safe fallback)
    output_list <- lapply(
      seq_len(mixtures),
      function(i) coreJob(i, Y = Y, X = X, absolute = absolute, 
                          abs_method = abs_method, Ymedian = Ymedian, 
                          P = P, nulldist = nulldist)
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
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE", paste("Absolute score (", abs_method, ")", sep = ""))
  }
  obj
}