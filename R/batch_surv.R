#' Batch Survival Analysis
#'
#' @description
#' Performs Cox proportional hazards regression analysis on multiple variables.
#' Optionally determines optimal cutoffs to dichotomize continuous predictors before
#' modeling. Returns hazard ratios, confidence intervals, and p-values for each variable.
#'
#' @param pdata Data frame containing survival time, event status, and predictor variables.
#' @param variable Character vector specifying the names of predictor variables to analyze.
#' @param time Character string specifying the column name containing follow-up time.
#'   Default is `"time"`.
#' @param status Character string specifying the column name containing event status
#'   (1 = event occurred, 0 = censored). Default is `"status"`.
#' @param best_cutoff Logical indicating whether to compute optimal cutoffs for
#'   continuous variables and analyze dichotomized versions. Default is `FALSE`.
#'
#' @return Data frame containing hazard ratios (HR), 95% confidence intervals (CI),
#'   and p-values for each variable, sorted by p-value.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' sig_stad <- load_data("sig_stad")
#' batch_surv(
#'   pdata = sig_stad,
#'   variable = colnames(sig_stad)[69:ncol(sig_stad)],
#'   time = "OS_time",
#'   status = "OS_status"
#' )
#' }
batch_surv <- function(pdata, variable, time = "time", status = "status", best_cutoff = FALSE) {
  if (!is.data.frame(pdata)) {
    cli::cli_abort("{.arg pdata} must be a data frame")
  }
  if (nrow(pdata) == 0) {
    cli::cli_abort("{.arg pdata} has no rows")
  }
  if (!time %in% colnames(pdata)) {
    cli::cli_abort("Time column {.val {time}} not found in pdata")
  }
  if (!status %in% colnames(pdata)) {
    cli::cli_abort("Status column {.val {status}} not found in pdata")
  }
  if (!is.character(variable) || length(variable) == 0) {
    cli::cli_abort("{.arg variable} must be a non-empty character vector")
  }

  pdata <- as.data.frame(pdata)
  colnames(pdata)[colnames(pdata) == time] <- "time_iobr"
  colnames(pdata)[colnames(pdata) == status] <- "status_iobr"

  pdata <- pdata[!is.na(pdata$time_iobr) & !is.na(pdata$status_iobr), , drop = FALSE]
  if (nrow(pdata) == 0) {
    cli::cli_abort("No valid data after removing NA values")
  }

  pdata$time_iobr <- as.numeric(pdata$time_iobr)
  pdata$status_iobr <- as.numeric(pdata$status_iobr)

  valid_vars <- variable[variable %in% colnames(pdata)]
  invalid_vars <- setdiff(variable, colnames(pdata))
  if (length(invalid_vars) > 0) {
    cli::cli_alert_warning("Ignoring {length(invalid_vars)} variable{?s} not found: {.val {invalid_vars}}")
  }
  if (length(valid_vars) == 0) {
    cli::cli_abort("No valid variables found in pdata")
  }

  if (best_cutoff) {
    cli::cli_alert_info("Computing optimal cutoffs for {length(valid_vars)} variable{?s}")
    for (var in valid_vars) {
      pdata <- best_cutoff(
        pdata = pdata,
        time = "time_iobr",
        status = "status_iobr",
        variable = var,
        print_result = FALSE
      )
      binary_col <- paste0(var, "_binary")
      pdata[[binary_col]] <- ifelse(pdata[[binary_col]] == "High", 1, 0)
    }
    valid_vars <- paste0(valid_vars, "_binary")
  }

  result_list <- lapply(valid_vars, function(var) {
    fit <- survival::coxph(
      survival::Surv(pdata$time_iobr, pdata$status_iobr) ~ pdata[[var]],
      data = pdata
    )
    getHRandCIfromCoxph(fit)
  })

  result <- dplyr::bind_rows(result_list)
  result$ID <- if (best_cutoff) gsub("_binary$", "", valid_vars) else valid_vars
  # 把 ID 调整到第一列
  result <- result[, c("ID", setdiff(colnames(result), "ID")), drop = FALSE]
  result <- result[order(result$P, decreasing = FALSE), , drop = FALSE]
  rownames(result) <- NULL

  tibble::as_tibble(result)
}
