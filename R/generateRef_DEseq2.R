#' Generate Reference Signature Matrix Using DESeq2
#'
#' This function uses DESeq2 to perform differential expression analysis across cell types,
#' identifies significantly expressed genes, and creates a reference signature matrix from median expression levels.
#'
#' @param dds Raw count data from RNA-seq (matrix).
#' @param pheno Character vector of cell type classes for samples.
#' @param FDR Numeric threshold for adjusted p-values to consider genes significant. Default is 0.05.
#' @param dat Normalized expression data (e.g., FPKM, TPM) for calculating median expression.
#' @return A list containing:
#'   - `reference_matrix`: Data frame of median expression for significant genes across cell types.
#'   - `G`: Optimal number of probes minimizing condition number.
#'   - `condition_number`: Minimum condition number.
#'   - `whole_matrix`: Full median expression matrix.
#' @export
#' @examples
#' # Example assuming 'dds' is raw counts and 'dat' is normalized data
#' dds <- matrix(sample(0:1000, 2000, replace = TRUE), nrow = 100, ncol = 20)
#' colnames(dds) <- paste("Sample", 1:20, sep = "_")
#' rownames(dds) <- paste("Gene", 1:100, sep = "_")
#' pheno <- rep(c("Type1", "Type2"), each = 10)
#' dat <- matrix(runif(2000), nrow = 100, ncol = 20)
#' rownames(dat) <- rownames(dds)
#' colnames(dat) <- colnames(dds)
#' result <- generateRef_DEseq2(dds = dds, pheno = pheno, FDR = 0.05, dat = dat)
#' print(result$reference_matrix)
generateRef_DEseq2 <- function(dds, pheno, FDR = 0.05, dat) {
  if (!all(colnames(dds) == colnames(dat))) {
    stop("The sample order of dds and dat much be identical")
  }
  dds <- round(dds, 0)
  keep <- rowSums(dds) >= 0.05 * ncol(dds)
  dds <- dds[keep, ]
  dds <- na.omit(dds)
  samples <- data.frame(row.names = colnames(dds), cell = pheno)
  ncells <- unique(pheno)
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
    purrr::map(function(x) Top_probe(dat = x, i = i)) %>%
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


#' Top Probe Selector
#'
#' This function extracts the top `i` probes based on their ordering in the provided data frame.
#' If the number of rows in the data frame is less than or equal to `i`, it returns all probes.
#'
#' @param dat A data frame containing a column named "probe" among other data.
#' @param i An integer indicating the number of top probes to return from the dataset.
#'
#' @return A character vector containing the names of the top `i` probes.
#' @export
#'
#' @examples
#' # Assuming 'dat' is a data frame with at least one column named "probe"
#' dat <- data.frame(
#'   probe = c("Probe1", "Probe2", "Probe3", "Probe4", "Probe5"),
#'   value = c(5, 3, 2, 4, 1)
#' )
#' # Get the top 3 probes
#' top_probes <- Top_probe(dat, 3)
#' print(top_probes)
Top_probe <- function(dat, i) {
  if (nrow(dat) <= i) {
    probe <- dat[, "probe"] %>% as.character()
  } else {
    probe <- dat[1:i, "probe"] %>% as.character()
  }
  probe
}
