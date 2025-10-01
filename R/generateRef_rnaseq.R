#' Generate Reference Gene Matrix from RNA-seq DEGs
#'
#' Uses DESeq2 to identify differentially expressed genes and create a reference matrix
#' from median expression levels across cell types.
#'
#' @param dds Raw count data from RNA-seq.
#' @param pheno Character vector of cell type classes.
#' @param mode "oneVSothers" or "pairs" for DEG identification.
#' @param FDR Numeric threshold for adjusted p-values. Default is 0.05.
#' @param dat Normalized expression data (e.g., FPKM, TPM).
#' @return A list containing reference matrix, optimal G, condition number, and whole matrix.
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @export
#' @author Rongfang Shen
#' @examples
#' dds <- matrix(rpois(200 * 10, lambda = 10), ncol = 10)
#' pheno <- sample(c("Type1", "Type2", "Type3"), 10, replace = TRUE)
#' dat <- matrix(rnorm(200 * 10), ncol = 10)
#' results <- generateRef_rnaseq(dds = dds, pheno = pheno, FDR = 0.05, dat = dat)
generateRef_rnaseq <- function(dds, pheno, mode = "oneVSothers", FDR = 0.05, dat) {
  if (!all(colnames(dds) == colnames(dat))) {
    stop("The sample order of dds and dat much be identical")
  }
  dds <- round(dds, 0)
  keep <- rowSums(dds) >= 0.05 * ncol(dds)
  dds <- dds[keep, ]
  dds <- na.omit(dds)
  samples <- data.frame(row.names = colnames(dds), cell = pheno)
  ncells <- unique(pheno)
  if (mode == "pairs") {
    dds <- DESeq2::DESeqDataSetFromMatrix(dds,
      colData = samples,
      design = ~cell
    )
    dds <- DESeq2::DESeq(dds)
    contrasts <- list()
    # generate contrast list
    for (i in 1:length(ncells)) {
      cellA <- ncells[i]
      othercells <- setdiff(ncells, cellA)
      z1 <- length(othercells)
      for (j in 1:length(othercells)) {
        cellB <- othercells[j]
        z <- (i - 1) * z1 + j
        contrasts[[z]] <- c("cell", cellA, cellB)
      }
    }
    res <- lapply(contrasts, function(z) {
      DESeq2::results(dds, contrast = z)
    })
  }
  if (mode == "oneVSothers") {
    res <- list()
    for (i in 1:length(ncells)) {
      cellA <- ncells[i]
      samples$group <- NA
      samples$group <- ifelse(samples$cell == cellA, cellA, "others")
      dds2 <- DESeq2::DESeqDataSetFromMatrix(dds,
        colData = samples,
        design = ~group
      )
      dds2 <- DESeq2::DESeq(dds2)
      res[[i]] <- DESeq2::results(dds2, contrast = c("group", cellA, "others"))
    }
  }
  resData <- lapply(res, function(x) {
    tmp <- as.data.frame(x[order(x$padj), ])
    tmp <- data.frame(probe = rownames(tmp), tmp)
    tmp <- tmp[tmp$padj < FDR, ]
    return(tmp)
  })
  median_value <- t(dat) %>%
    as.data.frame() %>%
    split(., pheno) %>%
    purrr::map(function(x) matrixStats::colMedians(as.matrix(x)))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- c()
  for (i in 50:200) {
    probes <- resData %>%
      purrr::map(function(x) Top_probe(dat = x, i = i)) %>%
      unlist() %>%
      unique()
    tmpdat <- median_value[probes, ]
    # condition number
    con_num <- kappa(tmpdat)
    con_nums <- c(con_nums, con_num)
  }
  i <- c(50:200)[which.min(con_nums)]
  probes <- resData %>%
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
