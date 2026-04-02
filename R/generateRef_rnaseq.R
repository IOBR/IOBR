#' Generate Reference Gene Matrix from RNA-seq DEGs
#'
#' @description
#' Uses DESeq2 to identify differentially expressed genes and create a
#' reference matrix from median expression levels across cell types.
#'
#' @param dds Matrix. Raw count data from RNA-seq.
#' @param pheno Character vector. Cell type classes.
#' @param mode Character. DEG identification mode: `"oneVSothers"` or
#'   `"pairs"`. Default is `"oneVSothers"`.
#' @param FDR Numeric. Threshold for adjusted p-values. Default is 0.05.
#' @param dat Matrix. Normalized expression data (e.g., FPKM, TPM).
#'
#' @return List containing:
#'   - `reference_matrix`: Data frame of median expression for significant
#'     genes across cell types.
#'   - `G`: Optimal number of probes minimizing condition number.
#'   - `condition_number`: Minimum condition number.
#'   - `whole_matrix`: Full median expression matrix.
#'
#' @export
#' @author Rongfang Shen
#'
#' @examples
#' dds <- matrix(rpois(200 * 10, lambda = 10), ncol = 10)
#' rownames(dds) <- paste("Gene", 1:200, sep = "_")
#' colnames(dds) <- paste("Sample", 1:10, sep = "_")
#' pheno <- sample(c("Type1", "Type2", "Type3"), 10, replace = TRUE)
#' dat <- matrix(rnorm(200 * 10), ncol = 10)
#' rownames(dat) <- rownames(dds)
#' colnames(dat) <- colnames(dds)
#' \donttest{
#' results <- generateRef_rnaseq(dds = dds, pheno = pheno, FDR = 0.05, dat = dat)
#' }
generateRef_rnaseq <- function(dds, pheno, mode = "oneVSothers", FDR = 0.05, dat) {
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

  mode <- rlang::arg_match(mode, c("oneVSothers", "pairs"))

  dds <- round(dds, 0)
  keep <- rowSums(dds) >= 0.05 * ncol(dds)
  dds <- dds[keep, ]
  dds <- stats::na.omit(dds)

  samples <- data.frame(row.names = colnames(dds), cell = pheno)
  ncells <- unique(pheno)

  if (mode == "pairs") {
    dds_obj <- DESeq2::DESeqDataSetFromMatrix(
      dds,
      colData = samples,
      design = ~cell
    )
    dds_obj <- DESeq2::DESeq(dds_obj)

    contrasts <- list()
    for (i in seq_along(ncells)) {
      cellA <- ncells[i]
      othercells <- setdiff(ncells, cellA)
      for (j in seq_along(othercells)) {
        cellB <- othercells[j]
        z <- (i - 1) * length(othercells) + j
        contrasts[[z]] <- c("cell", cellA, cellB)
      }
    }
    res <- lapply(contrasts, function(z) DESeq2::results(dds_obj, contrast = z))
  }

  if (mode == "oneVSothers") {
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
  }

  resData <- lapply(res, function(x) {
    tmp <- as.data.frame(x[order(x$padj), ])
    tmp <- data.frame(probe = rownames(tmp), tmp)
    tmp <- tmp[tmp$padj < FDR, ]
    tmp
  })

  median_value <- t(dat) |>
    as.data.frame() |>
    split(pheno) |>
    purrr::map(function(x) apply(as.matrix(x), 2, stats::median, na.rm = TRUE))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- numeric(151)
  for (i in 50:200) {
    probes <- resData |>
      purrr::map(function(x) Top_probe(dat = x, i = i)) |>
      unlist() |>
      unique()
    tmpdat <- median_value[probes, , drop = FALSE]
    con_nums[i - 49] <- kappa(tmpdat)
  }

  i <- 49 + which.min(con_nums)
  probes <- resData |>
    purrr::map(function(x) Top_probe(dat = x, i = i)) |>
    unlist() |>
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
