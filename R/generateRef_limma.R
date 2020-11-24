#'  Generate Reference signature matrix using limma
#'
#' @param dat data frame or matrix; gene probes in the row and samples in the column
#' @param pheno character vector; cell type class of the samples
#' @param FDR numeric; genes with BH adjust p value < FDR are considered significant.
#' @import tidyverse
#' @import stringr
#' @import purrr
#'
#' @return
#' @export
#'
#' @examples
generateRef_limma <- function(dat, pheno, FDR = 0.05){
  pheno <- factor(pheno)
  design <- model.matrix(~0 + pheno)
  colnames(design) <- stringr:: str_sub(colnames(design), 6)
  rownames(design) <- colnames(dat)
  con <- Construct_con(pheno = pheno)
  fit <- limma::lmFit(dat, design) %>% limma::contrasts.fit(., con) %>% limma::eBayes()
  dif <- list()
  dif <- colnames(con) %>%
    purrr::map(function(x) limma::topTable(fit, adjust.method="BH", coef = x, number = Inf)) %>%
    lapply(., function(x)data.frame(probe = rownames(x), x)) %>%
    purrr::map(function(x)filter(x,adj.P.Val < FDR)) %>%
    purrr::map(function(x)arrange(x, desc(logFC)))

  median_value <- t(dat) %>% as.data.frame() %>%
    split(., pheno) %>% map(function(x) matrixStats:: colMedians(as.matrix(x)))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- c()
  for (i in 50:200){
    probes <- dif %>% map(function(x)Top_probe(dat = x, i = i)) %>%
      unlist() %>% unique()
    tmpdat <- median_value[probes, ]
    #condition number
    con_num <- kappa(tmpdat)
    con_nums <- c(con_nums, con_num)
  }
  i <- c(50:200)[which.min(con_nums)]
  probes <- dif %>% map(function(x)Top_probe(dat = x, i = i)) %>%
    unlist() %>% unique()
  reference <- median_value[probes, ]
  reference <- data.frame(NAME = rownames(reference), reference)
  G <- i
  condition_number <- min(con_nums)

  return(list(reference_matrix = reference, G = G, condition_number = condition_number,
              whole_matrix = median_value))
}



#' Construct contrast
#'
#' @param pheno
#' @param mode
#'
#' @return
#' @export
#'
#' @examples
Construct_con <- function(pheno, mode){
  n <- length(levels(pheno))
  con <- matrix(-1, n, n)
  diag(con) <- 1
  rownames(con) <- levels(pheno)
  colnames(con) <- paste0(rownames(con), " - others")
  return(con)
}
