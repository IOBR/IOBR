#' Generate Reference Signature Matrix Using Limma
#'
#' This function performs differential expression analysis using the limma package to identify significantly
#' expressed genes across different cell types specified in `pheno`. It computes median expression levels
#' of these significant genes to create a reference signature matrix. The function is particularly useful
#' for constructing signature matrices in gene expression studies involving multiple cell types or conditions.
#'
#'
#' @param dat data frame or matrix; gene probes in the row and samples in the column
#' @param pheno character vector; cell type class of the samples
#' @param FDR numeric; genes with BH adjust p value < FDR are considered significant.
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import stringr
#' @import purrr
#'
#' @return A list containing:
#'         - `reference_matrix`: A dataframe of median expression values for significantly expressed genes, with genes as rows.
#'         - `G`: The number of probes used in the reference matrix that resulted in the lowest condition number.
#'         - `condition_number`: The minimum condition number obtained from the analysis.
#'         - `whole_matrix`: A matrix of median values across all samples and conditions for further analysis.
#' @export
#'
#' @examples
#' # Assume 'dat' is a matrix of gene expression data and 'pheno' is the corresponding cell type information
#' dat <- matrix(rnorm(2000), nrow = 100)
#' pheno <- sample(c("Type1", "Type2", "Type3"), 20, replace = TRUE)
#' results <- generateRef_limma(dat, pheno)
#' print(results)
generateRef_limma <- function(dat, pheno, FDR = 0.05) {
  pheno <- factor(pheno)
  design <- model.matrix(~ 0 + pheno)
  colnames(design) <- stringr::str_sub(colnames(design), 6)
  rownames(design) <- colnames(dat)
  con <- Construct_con(pheno = pheno)
  fit <- limma::lmFit(dat, design) %>%
    limma::contrasts.fit(., con) %>%
    limma::eBayes()
  dif <- list()
  dif <- colnames(con) %>%
    purrr::map(function(x) limma::topTable(fit, adjust.method = "BH", coef = x, number = Inf)) %>%
    lapply(., function(x) data.frame(probe = rownames(x), x)) %>%
    purrr::map(function(x) filter(x, adj.P.Val < FDR)) %>%
    purrr::map(function(x) arrange(x, desc(logFC)))

  median_value <- t(dat) %>%
    as.data.frame() %>%
    split(., pheno) %>%
    map(function(x) matrixStats::colMedians(as.matrix(x)))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- c()
  for (i in 50:200) {
    probes <- dif %>%
      map(function(x) Top_probe(dat = x, i = i)) %>%
      unlist() %>%
      unique()
    tmpdat <- median_value[probes, ]
    # condition number
    con_num <- kappa(tmpdat)
    con_nums <- c(con_nums, con_num)
  }
  i <- c(50:200)[which.min(con_nums)]
  probes <- dif %>%
    map(function(x) Top_probe(dat = x, i = i)) %>%
    unlist() %>%
    unique()
  reference <- median_value[probes, ]
  reference <- data.frame(NAME = rownames(reference), reference)
  G <- i
  condition_number <- min(con_nums)

  return(list(
    reference_matrix = reference, G = G, condition_number = condition_number,
    whole_matrix = median_value
  ))
}



#' Construct Contrast Matrix
#'
#' This function creates a contrast matrix for differential analysis, where each phenotype level is contrasted
#' against all other levels combined. This is particularly useful in linear models for comparing multiple groups.
#'
#' @param pheno A factor variable with different levels representing groups or conditions to contrast.
#' @param mode Currently unused but reserved for future extensions where different modes of contrast might be implemented.
#'
#' @return A square matrix with dimensions equal to the number of levels in `pheno`.
#'         Each row represents a contrast where the corresponding level is compared against the average of others.
#'         The matrix elements are set to -1 for non-diagonal cells (indicating comparison groups)
#'         and 1 for diagonal cells (indicating the group of interest).
#' @export
#'
#' @examples
#' # Example usage:
#' pheno <- factor(c("A", "B", "C", "D"))
#' contrast_matrix <- Construct_con(pheno, mode = NULL)
#' print(contrast_matrix)
Construct_con <- function(pheno, mode) {
  n <- length(levels(pheno))
  con <- matrix(-1, n, n)
  diag(con) <- 1
  rownames(con) <- levels(pheno)
  colnames(con) <- paste0(rownames(con), " - others")
  return(con)
}
