#' Generate reference gene matrix
#'
#' @param dat data frame or matrix; gene probes in the row and samples in the column
#' @param pheno character vector; cell type class of the samples
#' @param mode "oneVSothers" or "pairs"; two modes for identifying significantly differential expressed genes
#' @param FDR numeric; genes with BH adjust p value < FDR are considered significant.
#'
#' @return
#' @export
#' @author Rongfang Shen
#'
#' @examples

generateRef <- function(dat, pheno, mode = c("oneVSothers", "pairs"), FDR = 0.05){
  pheno <- factor(pheno)
  design <- model.matrix(~0 + pheno)
  colnames(design) <- stringr:: str_sub(colnames(design), 6)
  rownames(design) <- colnames(dat)
  con <- Construct_con(pheno = pheno, mode = mode)
  fit <- lmFit(dat, design) %>% contrasts.fit(., con) %>% eBayes()
  dif <- list()
  dif <- colnames(con) %>%
    purrr::map(function(x) limma::topTable(fit, adjust.method="BH", coef = x, number = Inf)) %>%
    lapply(., function(x) data.frame(probe = rownames(x), x)) %>%
    purrr::map(function(x)dplyr:: filter(x, adj.P.Val < FDR)) %>%
    purrr::map(function(x)dplyr:: arrange(x, desc(logFC)))

  median_value <- t(dat) %>% as.data.frame() %>%
    split(., pheno) %>% purrr:: map(function(x)matrixStats::colMedians(as.matrix(x)))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  # select top G genes: 50-200; condition number
  con_nums <- c()
  for (i in 50:200){
    probes <- dif %>% purrr::map(function(x) Top_probe(dat = x, i = i)) %>%
      unlist() %>% unique()
    tmpdat <- median_value[probes, ]
    #condition number
    con_num <- kappa(tmpdat)
    con_nums <- c(con_nums, con_num)
  }
  i <- c(50:200)[which.min(con_nums)]
  probes <- dif %>% purrr::map(function(x)Top_probe(dat = x, i = i)) %>%
    unlist() %>% unique()
  reference <- median_value[probes, ]
  reference <- data.frame(NAME = rownames(reference), reference)
  G <- i
  condition_number <- min(con_nums)

  return(list(reference_matrix = reference, G = G, condition_number = condition_number,
              whole_matrix = median_value))
}


#' Title
#'
#' @param dat
#' @param i
#'
#' @return
#' @export
#'
#' @examples
Top_probe <- function(dat, i){
  if (nrow(dat) <= i){
    probe = dat[, "probe"] %>% as.character()
  }else{
    probe = dat[1:i, "probe"] %>% as.character()
  }
  probe
}


#' Construct contrast formula
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
  if (mode == "oneVSothers"){
    con <- matrix(-1, n, n)
    diag(con) <- 1
    rownames(con) <- levels(pheno)
    colnames(con) <- paste0(rownames(con), " - others")
  }
  if (mode == "pairs"){
    con <- matrix(0, n, choose(n,2))
    rownames(con) <- levels(pheno)
    colnames(con) <- 1:choose(n,2)
    k = 0
    for (i in 1:(n -1)){
      for (j in (i+1):n){
        k = k + 1
        con[i, k] = 1
        con[j, k] = -1
        colnames(con)[k] <- paste0(rownames(con)[i], "-", rownames(con)[j])
      }
    }
  }
  return(con)
}
