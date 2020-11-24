#' Generate Reference signature matrix using DEseq2
#'
#' @param dds raw count data from RNA-seq
#' @param pheno  character vector; cell type class of the samples
#' @param FDR numeric; genes with BH adjust p value < FDR are considered significant.
#' @param dat data frame or matrix; normalized transcript quantification data (like FPKM, TPM). Note: cell's median expression level of the identified probes will be the output of reference_matrix.
#' @return
#' @export
#' @example
generateRef_DEseq2 <- function(dds, pheno, FDR = 0.05, dat){
  if (!all(colnames(dds) == colnames(dat))){
    stop("The sample order of dds and dat much be identical")
  }
  dds <- round(dds, 0)
  keep <- rowSums(dds) >= 0.05 *ncol(dds)
  dds <- dds[keep,]
  dds <- na.omit(dds)
  samples <- data.frame(row.names = colnames(dds), cell = pheno)
  ncells <- unique(pheno)
  res <- list()
  for (i in 1:length(ncells)){
     cellA <- ncells[i]
     samples$group <- NA
     samples$group <- ifelse(samples$cell == cellA, cellA, "others")
     dds2 <- DESeq2::DESeqDataSetFromMatrix(dds,
                                   colData = samples,
                                   design = ~ group)
     dds2 <- DESeq2::DESeq(dds2)
     res[[i]] <- DESeq2::results(dds2, contrast = c("group", cellA, "others"))
    }

  resData <- lapply(res, function(x){
    tmp <- as.data.frame(x[order(x$padj), ])
    tmp <- data.frame(probe = rownames(tmp), tmp)
    tmp <- tmp[tmp$padj < FDR, ]
    return(tmp)
  })
  median_value <- t(dat) %>% as.data.frame() %>%
    split(., pheno) %>% purrr::map(function(x) matrixStats:: colMedians(as.matrix(x)))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- c()
  for (i in 50:200){
    probes <- resData %>% purrr::map(function(x)Top_probe(dat = x, i = i)) %>%
      unlist() %>% unique()
    tmpdat <- median_value[probes, ]
    #condition number
    con_num <- kappa(tmpdat)
    con_nums <- c(con_nums, con_num)
  }
  i <- c(50:200)[which.min(con_nums)]
  probes <- resData %>% purrr::map(function(x)Top_probe(dat = x, i = i)) %>%
    unlist() %>% unique()
  reference <- median_value[probes, ]
  reference <- data.frame(NAME = rownames(reference), reference)
  G <- i
  condition_number <- min(con_nums)
  return(list(reference_matrix = reference, G = G, condition_number = condition_number,
              whole_matrix = median_value))
}


#' Top probe
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
