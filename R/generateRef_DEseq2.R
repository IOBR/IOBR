#' Generate Reference Signature Matrix Using DESeq2
#'
#' @description
#' Uses DESeq2 to perform differential expression analysis across cell types,
#' identifies significantly expressed genes, and creates a reference signature
#' matrix from median expression levels.
#'
#' @param dds Matrix. Raw count data from RNA-seq.
#' @param pheno Character vector. Cell type classes for samples.
#' @param FDR Numeric. Threshold for adjusted p-values. Default is 0.05.
#' @param dat Matrix. Normalized expression data (e.g., FPKM, TPM) for
#'   calculating median expression.
#'
#' @return List containing:
#'   - `reference_matrix`: Data frame of median expression for significant
#'     genes across cell types.
#'   - `G`: Optimal number of probes minimizing condition number.
#'   - `condition_number`: Minimum condition number.
#'   - `whole_matrix`: Full median expression matrix.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' dds <- matrix(sample(0:1000, 2000, replace = TRUE), nrow = 100, ncol = 20)
#' colnames(dds) <- paste("Sample", 1:20, sep = "_")
#' rownames(dds) <- paste("Gene", 1:100, sep = "_")
#' pheno <- rep(c("Type1", "Type2"), each = 10)
#' dat <- matrix(runif(2000), nrow = 100, ncol = 20)
#' rownames(dat) <- rownames(dds)
#' colnames(dat) <- colnames(dds)
#' \donttest{
#' result <- generateRef_DEseq2(dds = dds, pheno = pheno, FDR = 0.05, dat = dat)
#' print(result$reference_matrix)
#' }
generateRef_DEseq2 <- function(dds, pheno, FDR = 0.05, dat) {
  rlang::check_installed("DESeq2")

  if (!is.matrix(dds) && !is.data.frame(dds)) {
    cli::cli_abort("{.arg dds} must be a matrix or data frame")
  }
  if (!all(colnames(dds) == colnames(dat))) {
    cli::cli_abort("Sample order of {.arg dds} and {.arg dat} must be identical")
  }
  if (length(pheno) != ncol(dds)) {
    cli::cli_abort("Length of {.arg pheno} must match number of columns in {.arg dds}")
  }

  dds <- round(dds, 0)
  keep <- rowSums(dds) >= 0.05 * ncol(dds)
  dds <- dds[keep, ]
  dds <- stats::na.omit(dds)

  samples <- data.frame(row.names = colnames(dds), cell = pheno)
  ncells <- unique(pheno)
  res <- vector("list", length(ncells))

  for (i in seq_along(ncells)) {
    cellA <- ncells[i]
    samples$group <- ifelse(samples$cell == cellA, cellA, "others")
    dds2 <- DESeq2::DESeqDataSetFromMatrix(
      dds,
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
    tmp
  })

  median_value <- t(dat) %>%
    as.data.frame() %>%
    split(pheno) %>%
    purrr::map(function(x) apply(as.matrix(x), 2, stats::median, na.rm = TRUE))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- numeric(151)
  for (i in 50:200) {
    probes <- resData %>%
      purrr::map(function(x) Top_probe(dat = x, i = i)) %>%
      unlist() %>%
      unique()
    tmpdat <- median_value[probes, , drop = FALSE]
    con_nums[i - 49] <- kappa(tmpdat)
  }

  i <- 49 + which.min(con_nums)
  probes <- resData %>%
    purrr::map(function(x) Top_probe(dat = x, i = i)) %>%
    unlist() %>%
    unique()
  reference <- median_value[probes, , drop = FALSE]
  reference <- data.frame(NAME = rownames(reference), reference)
  G <- i
  condition_number <- min(con_nums)

  list(
    reference_matrix = reference,
    G = G,
    condition_number = condition_number,
    whole_matrix = median_value
  )
}


#' Top Probe Selector
#'
#' @description
#' Extracts the top `i` probes based on their ordering in the provided data
#' frame. If the number of rows is less than or equal to `i`, returns all
#' probes.
#'
#' @param dat Data frame containing a column named "probe".
#' @param i Integer. Number of top probes to return.
#'
#' @return Character vector containing the names of the top `i` probes.
#'
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   probe = c("Probe1", "Probe2", "Probe3", "Probe4", "Probe5"),
#'   value = c(5, 3, 2, 4, 1)
#' )
#' top_probes <- Top_probe(dat, 3)
#' print(top_probes)
Top_probe <- function(dat, i) {
  if (nrow(dat) <= i) {
    probe <- as.character(dat[, "probe"])
  } else {
    probe <- as.character(dat[seq_len(i), "probe"])
  }
  probe
}
