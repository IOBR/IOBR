#' Title
#'
#' @param dat data frame or matrix; gene probes in the row and samples in the column
#' @param pheno character vector; cell type class of the samples
#' @param mode "oneVSothers" or "pairs"; two modes for identifying significantly differentially expressed genes
#' @param FDR numeric; genes with BH adjust p value < FDR are considered significant.
#' @import limma
#' @import tidyverse
#' @import matrixStats
#' @import purrr
#'
#' @return
#' @export
#'
#' @examples
#' data("Stromal_v2.pheno")
#' data("Stromal_v2")
#' immune_feature <- generateRef(dat = Stromal_v2, pheno = Stromal_v2.pheno, mode = "oneVSothers", FDR = 0.05)
#' reference <- immune_feature$reference_matrix
#' condition_number <- immune_feature$condition_number
#' top_probe <- immune_feature$G
generateRef <- function(dat, pheno, mode = c("oneVSothers", "pairs"), FDR){
  pheno <- factor(pheno)
  design <- model.matrix(~0 + pheno)
  colnames(design) <- str_sub(colnames(design), 6)
  rownames(design) <- colnames(dat)
  con <- Construct_con(pheno = pheno, mode = mode)
  fit <- lmFit(dat, design) %>% contrasts.fit(., con) %>% eBayes()
  dif <- list()
  dif <- colnames(con) %>%
    purrr::map(function(x) limma::topTable(fit, adjust.method="BH", coef = x, number = Inf)) %>%
    lapply(., function(x)data.frame(probe = rownames(x), x)) %>%
    purrr::map(function(x)filter(x,adj.P.Val < FDR)) %>%
    purrr::map(function(x)arrange(x, desc(logFC)))

  median_value <- t(dat) %>% as.data.frame() %>%
    split(., pheno) %>% map(function(x)colMedians(as.matrix(x)))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  # select top G genes: 50-200; condition number
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

Top_probe <- function(dat, i){
  if (nrow(dat) <= i){
    probe = dat[, "probe"] %>% as.character()
  }else{
    probe = dat[1:i, "probe"] %>% as.character()
  }
  probe
}
Construct_con <- function(pheno, mode){
  n <- length(levels(pheno))
  if (mode == "oneVSothers"){
    con <- matrix(-1, n, n)
    diag(con) <- 1
    rownames(con) <- levels(pheno)
    colnames(con) <- paste0(rownames(con), " - others")
  }
  if (mode == "pairs"){
    con <- matrix(0,n,choose(n,2))
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
