#' Combine Phenotype Data and Expression Set
#'
#' @description
#' Merges phenotype data with an expression matrix by matching sample IDs.
#' Optionally filters features, applies feature manipulation, and scales
#' expression data before combining.
#'
#' @param eset Expression matrix with genes/features in rows and samples in
#'   columns.
#' @param pdata Data frame containing phenotype/clinical data.
#' @param id_pdata Character string specifying the column name in `pdata`
#'   containing sample identifiers. Default is `"ID"`.
#' @param feas Character vector specifying features to include from `eset`.
#'   If `NULL`, all features are used. Default is `NULL`.
#' @param feature_manipulation Logical indicating whether to apply feature
#'   manipulation to filter valid features. Default is `TRUE`.
#' @param scale Logical indicating whether to scale (standardize) expression
#'   data. Default is `TRUE`.
#' @param choose_who_when_duplicate Character string specifying which data to
#'   prefer when duplicate columns exist. Options are `"eset"` or `"pdata"`.
#'   Default is `"eset"`.
#'
#' @return Data frame combining phenotype data and (transposed) expression data,
#'   with samples in rows and features/phenotypes in columns.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' set.seed(123)
#' eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(eset) <- paste0("Gene", 1:100)
#' colnames(eset) <- paste0("Sample", 1:10)
#' pdata <- data.frame(
#'   ID = colnames(eset),
#'   group = rep(c("A", "B"), each = 5),
#'   age = rnorm(10, 50, 10)
#' )
#' result <- combine_pd_eset(eset = eset, pdata = pdata, scale = FALSE)
#' dim(result)
combine_pd_eset <- function(eset,
                            pdata,
                            id_pdata = "ID",
                            feas = NULL,
                            feature_manipulation = TRUE,
                            scale = TRUE,
                            choose_who_when_duplicate = c("eset", "pdata")) {
  # Validate arguments
  choose_who_when_duplicate <- rlang::arg_match(choose_who_when_duplicate)

  # Filter features if specified
  if (!is.null(feas)) {
    feas <- feas[feas %in% rownames(eset)]
    if (length(feas) == 0) {
      cli::cli_abort("No specified features found in eset")
    }
    eset <- eset[feas, , drop = FALSE]
    cli::cli_alert_info("Filtered to {nrow(eset)} features")
  }

  # Apply feature manipulation
  if (feature_manipulation) {
    feas_filter <- IOBR::feature_manipulation(
      data = eset,
      feature = rownames(eset),
      is_matrix = TRUE,
      print_result = FALSE
    )
    eset <- eset[feas_filter, , drop = FALSE]
    cli::cli_alert_info("Retained {nrow(eset)} features after QC")
  }

  # Transpose and scale
  eset <- t(eset)
  if (scale) {
    eset <- scale(eset, center = TRUE, scale = TRUE)
  }

  # Convert to data frame with ID column
  eset_df <- as.data.frame(eset)
  eset_df$ID <- rownames(eset_df)

  # Check pdata ID column
  if (!id_pdata %in% colnames(pdata)) {
    cli::cli_abort("Column {.val {id_pdata}} not found in pdata")
  }
  pdata_df <- as.data.frame(pdata)
  colnames(pdata_df)[colnames(pdata_df) == id_pdata] <- "ID"

  # Merge data frames
  choose <- ifelse(choose_who_when_duplicate == "eset", "x", "y")

  pd_eset <- merge_duplicate(
    pdata_df,
    eset_df,
    by.x = "ID",
    by.y = "ID",
    choose = choose,
    all = FALSE
  )

  cli::cli_alert_success(
    "Combined data: {nrow(pd_eset)} samples x {ncol(pd_eset)} variables"
  )

  pd_eset
}
