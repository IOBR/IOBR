#' Check Integrity and Outliers of Expression Set
#'
#' @description
#' Performs quality checks on an expression matrix to identify missing values, infinite
#' values, and features with zero variance. Issues warnings when potential problems are
#' detected that may affect downstream analyses.
#'
#' @param eset Expression matrix with genes/features in rows and samples in columns.
#' @param print_result Logical indicating whether to print detailed check results to
#'   the console. Default is \code{FALSE}.
#' @param estimate_sd Logical indicating whether to check for features with zero
#'   standard deviation. Default is \code{FALSE}.
#'
#' @return Invisibly returns \code{NULL}. Function is called for its side effects
#'   (printing messages and issuing warnings).
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD expression data
#' data("eset_stad", package = "IOBR")
#' # Convert counts to TPM
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' # Check expression set integrity
#' check_eset(eset)
check_eset <- function(eset, print_result = FALSE, estimate_sd = FALSE) {
  if (print_result) {
    message(paste("<<< Do NA values exist ? >>> "))
    print(sum(is.na(eset)))
  }

  if (sum(is.na(eset)) > 0) warning(paste0("There are some missing values in the matrix, which may affect the score calculation. You can set parameter 'adjust_eset' as TRUE to avoid these effects"))


  teset <- as.data.frame(t(eset))

  if (print_result) {
    message(paste0("<<< Do -Inf features exist ? >>>"))
    print(summary(lapply(teset, function(x) min(x)) == -Inf))
  }


  if (print_result) {
    message(paste0("<<< Do Inf features exist ? >>>"))
    print(summary(lapply(teset, function(x) max(x)) == Inf))
  }

  if (min(eset, na.rm = TRUE) == -Inf | max(eset, na.rm = TRUE) == Inf) warning(paste0("There are infinite values in the matrix, which may affect the score calculation. You can set the 'adjust_eset' parameter as TRUE to avoid these effects."))


  if (estimate_sd) {
    sd <- apply(eset, 1, function(x) sd(x) == 0)

    if (print_result) {
      message(paste0("<<< Features have sd = 0 >>> "))
      print(summary(sd))
    }

    if (nlevels(as.factor(sd)) > 1) warning(paste0("Some variables in the matrix have no variance between samples, which may affect the score calculation. You can set the 'adjust_eset' parameter as TRUE to avoid these effects."))
  }
}
