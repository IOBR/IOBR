#' Generate Reference Signature Matrix Using Limma
#'
#' @description
#' Performs differential expression analysis using the limma package to
#' identify significantly expressed genes across different cell types.
#' Computes median expression levels of these significant genes to create
#' a reference signature matrix.
#'
#' @param dat Matrix or data frame. Gene probes in rows and samples in columns.
#' @param pheno Character vector. Cell type class of the samples.
#' @param FDR Numeric. Genes with BH adjusted p-value < FDR are considered
#'   significant. Default is 0.05.
#'
#' @return List containing:
#'   - `reference_matrix`: Data frame of median expression values for
#'     significantly expressed genes.
#'   - `G`: Number of probes used that resulted in the lowest condition number.
#'   - `condition_number`: Minimum condition number obtained.
#'   - `whole_matrix`: Matrix of median values across all samples.
#'
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(2000), nrow = 100)
#' rownames(dat) <- paste("Gene", 1:100, sep = "_")
#' colnames(dat) <- paste("Sample", 1:20, sep = "_")
#' pheno <- sample(c("Type1", "Type2", "Type3"), 20, replace = TRUE)
#' \donttest{
#' results <- generateRef_limma(dat, pheno)
#' print(results)
#' }
generateRef_limma <- function(dat, pheno, FDR = 0.05) {
  rlang::check_installed("limma")

  if (!is.matrix(dat) && !is.data.frame(dat)) {
    cli::cli_abort("{.arg dat} must be a matrix or data frame")
  }
  if (length(pheno) != ncol(dat)) {
    cli::cli_abort("Length of {.arg pheno} must match number of columns in {.arg dat}")
  }

  pheno <- factor(pheno)
  design <- stats::model.matrix(~ 0 + pheno)
  colnames(design) <- stringr::str_sub(colnames(design), 6)
  rownames(design) <- colnames(dat)

  con <- Construct_con(pheno = pheno)

  fit <- limma::lmFit(dat, design) |>
    limma::contrasts.fit(con) |>
    limma::eBayes()

  dif <- colnames(con) |>
    purrr::map(function(x) limma::topTable(fit, adjust.method = "BH", coef = x, number = Inf)) |>
    lapply(function(x) data.frame(probe = rownames(x), x)) |>
    purrr::map(function(x) dplyr::filter(x, .data$adj.P.Val < FDR)) |>
    purrr::map(function(x) dplyr::arrange(x, desc(.data$logFC)))

  median_value <- t(dat) |>
    as.data.frame() |>
    split(pheno) |>
    purrr::map(function(x) apply(as.matrix(x), 2, stats::median, na.rm = TRUE))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- numeric(151)
  for (i in 50:200) {
    probes <- dif |>
      purrr::map(function(x) Top_probe(dat = x, i = i)) |>
      unlist() |>
      unique()
    tmpdat <- median_value[probes, , drop = FALSE]
    con_nums[i - 49] <- kappa(tmpdat)
  }

  i <- 49 + which.min(con_nums)
  probes <- dif |>
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


#' Construct Contrast Matrix
#'
#' @description
#' Creates a contrast matrix for differential analysis, where each phenotype
#' level is contrasted against all other levels combined.
#'
#' @param pheno Factor with different levels representing groups to contrast.
#'
#' @return Square matrix with dimensions equal to the number of levels in
#'   `pheno`. Each row represents a contrast where the corresponding level
#'   is compared against the average of others.
#'
#' @export
#'
#' @examples
#' pheno <- factor(c("A", "B", "C", "D"))
#' contrast_matrix <- Construct_con(pheno)
#' print(contrast_matrix)
Construct_con <- function(pheno) {
  n <- length(levels(pheno))
  con <- matrix(-1, n, n)
  diag(con) <- 1
  rownames(con) <- levels(pheno)
  colnames(con) <- paste0(rownames(con), " - others")
  con
}
