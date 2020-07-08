#' Title
#'
#' @param dds raw count data from RNA-seq
#' @param pheno  data frame or matrix;un-normalized counts of RNA-seq
#' @param mode mode "oneVSothers" or "pairs"; two modes for identifying significantly differentially expressed genes
#' @param FDR character vector; cell type class of the samples
#' @param dat data frame or matrix; normalized transcript quantification data (like FPKM, TPM). Note: cell's median expression level of the identified probes will be the output of reference_matrix.
#' @import DESeq2
#'
#' @return
#' @export
#'
#' @example
generateRef_rnaseq <- function(dds, pheno, mode = "oneVSothers", FDR = 0.05, dat){
  if (!all(colnames(dds) == colnames(dat))){
    stop("The sample order of dds and dat much be identical")
  }
  dds <- round(dds, 0)
  keep <- rowSums(dds) >= 10
  dds <- dds[keep,]
  dds <- na.omit(dds)
  samples <- data.frame(row.names = colnames(dds), cell = pheno)
  ncells <- unique(pheno)
  if (mode == "pairs"){
    dds <- DESeqDataSetFromMatrix(dds,
                                    colData = samples,
                                    design = ~ cell)
    dds <- DESeq(dds)
    contrasts <- list()
    # generate contrast list
    for (i in 1:length(ncells)){
      cellA <- ncells[i]
      othercells <- setdiff(ncells, cellA)
      z1 <- length(othercells)
      for (j in 1:length(othercells)){
        cellB <- othercells[j]
        z <- (i - 1) * z1 + j
        contrasts[[z]]   <- c("cell", cellA, cellB)
      }
    }
    res <- lapply(contrasts, function(z){
      results(dds, contrast = z)
    })
  }
  if (mode == "oneVSothers"){
    res <- list()
    for (i in 1:length(ncells)){
     cellA <- ncells[i]
     samples$group <- NA
     samples$group <- ifelse(samples$cell == cellA, cellA, "others")
     dds2 <- DESeqDataSetFromMatrix(dds,
                                   colData = samples,
                                   design = ~ group)
     dds2 <- DESeq(dds2)
     res[[i]] <- results(dds2, contrast = c("group", cellA, "others"))
    }
  }
  resData <- lapply(res, function(x){
    tmp <- as.data.frame(x[order(x$padj), ])
    tmp <- data.frame(probe = rownames(tmp), tmp)
    tmp <- tmp[tmp$padj < FDR, ]
    return(tmp)
  })
  median_value <- t(dat) %>% as.data.frame() %>%
    split(., pheno) %>% map(function(x)colMedians(as.matrix(x)))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- c()
  for (i in 50:200){
    probes <- resData %>% map(function(x)Top_probe(dat = x, i = i)) %>%
      unlist() %>% unique()
    tmpdat <- median_value[probes, ]
    #condition number
    con_num <- kappa(tmpdat)
    con_nums <- c(con_nums, con_num)
  }
  i <- c(50:200)[which.min(con_nums)]
  probes <- resData %>% map(function(x)Top_probe(dat = x, i = i)) %>%
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
