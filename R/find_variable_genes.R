" Identify Variable Genes in Expression Data
#"
#' @description
#' Identifies variable genes from a gene expression dataset using specified
#' selection criteria. Supports multiple methods, including expression
#' thresholding and variability estimation via median absolute deviation (MAD).
#'
#' @param eset Numeric matrix. Gene expression data (genes as rows,
#'   samples as columns).
#' @param data_type Character. Type of data: "count" or "normalized".
#'   Default is "count".
#' @param methods Character vector. Methods for gene selection: "low", "mad".
#'   Default is c("low", "mad").
#' @param prop Numeric. Proportion of samples in which a gene must be expressed.
#'   Default is 0.7.
#' @param quantile Numeric. Quantile threshold for minimum MAD (0.25, 0.5, 0.75).
#'   Default is 0.75.
#' @param min.mad Numeric. Minimum allowable MAD value. Default is 0.1.
#' @param feas Character vector or NULL. Additional features to include.
#'   Default is NULL.
#'
#' @return Matrix subset of `eset` containing variable genes.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' eset_tme_stad <- load_data("eset_tme_stad")
#' eset <- find_variable_genes(
#'   eset = eset_tme_stad,
#'   data_type = "normalized",
#'   methods = "mad",
#'   quantile = 0.25
#' )
#' }
find_variable_genes <- function(eset,
                                data_type = c("count", "normalized"),
                                methods = c("low", "mad"),
                                prop = 0.7,
                                quantile = c(0.75, 0.5, 0.25),
                                min.mad = 0.1,
                                feas = NULL) {
  data_type <- rlang::arg_match(data_type)
  methods <- match.arg(methods, c("low", "mad"), several.ok = TRUE)
  quantile <- match.arg(as.character(quantile), c("0.25", "0.5", "0.75"))
  quantile <- as.numeric(quantile)

  eset <- as.matrix(eset)
  feas0 <- feature_manipulation(data = eset, is_matrix = TRUE)
  eset <- eset[feas0, , drop = FALSE]

  # Low expression filtering
  feas1 <- NULL
  if (data_type == "count" && "low" %in% methods) {
    cli::cli_alert_info(
      "Genes expressed in more than {prop * 100}% of samples"
    )
    keep <- rowSums(eset == 0) < ncol(eset) * prop
    print(table(keep))
    feas1 <- rownames(eset)[keep]
  }

  # MAD-based filtering
  feas2 <- NULL
  if ("mad" %in% methods) {
    eset_log <- log2eset(eset)
    cli::cli_alert_info("min.mad = {min.mad}")
    m.mad <- apply(eset_log, 1, stats::mad)
    cli::cli_alert_info("Range of MAD: {paste(round(range(m.mad), 2), collapse = ' to ')}")

    # Map quantile to index
    index <- switch(as.character(quantile),
      "0.75" = 4,
      "0.5" = 3,
      "0.25" = 2
    )

    cli::cli_alert_info("{quantile * 100}% of variables will be filtered out...")
    mad_quantiles <- quantile(m.mad, probs = seq(0, 1, 0.25))
    threshold <- max(mad_quantiles[index], min.mad)
    feas2 <- rownames(eset_log)[which(m.mad > threshold)]
  }

  feas <- unique(c(feas1, feas2, feas))
  feas <- feas[!is.na(feas)]

  eset[rownames(eset) %in% feas, , drop = FALSE]
}
