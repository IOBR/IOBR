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
#' # Load TCGA-STAD microenvironment signature data
#' sig_stad <- load_data("sig_stad")
#' # Perform batch survival analysis
#' batch_surv(
#'   pdata = sig_stad,
#'   variable = colnames(sig_stad)[69:ncol(sig_stad)],
#'   time = "OS_time",
#'   status = "OS_status"
#' )
#' }
batch_surv <- function(pdata, variable, time = "time", status = "status", best_cutoff = FALSE) {

  pdata <- as.data.frame(pdata)

  # Validate columns
  if (!time %in% colnames(pdata)) {
    cli::cli_abort("Time column {.val {time}} not found in pdata")
  }
  if (!status %in% colnames(pdata)) {
    cli::cli_abort("Status column {.val {status}} not found in pdata")
  }

  # Standardize column names
  colnames(pdata)[colnames(pdata) == time] <- "time_iobr"
  colnames(pdata)[colnames(pdata) == status] <- "status_iobr"

  # Remove NA values
  pdata <- pdata[!is.na(pdata$time_iobr) & !is.na(pdata$status_iobr), ]

  if (nrow(pdata) == 0) {
    cli::cli_abort("No valid data after removing NA values")
  }

  pdata$time_iobr <- as.numeric(pdata$time_iobr)
  pdata$status_iobr <- as.numeric(pdata$status_iobr)

  # Apply best cutoff if requested
  if (best_cutoff) {
    for (i in seq_along(variable)) {
      var <- variable[i]
      cli::cli_alert_info("Processing variable: {.val {var}}")

      pdata <- best_cutoff(
        pdata = pdata,
        time = "time_iobr",
        status = "status_iobr",
        variable = var,
        PrintResult = FALSE
      )

      binary_col <- paste0(var, "_binary")
      pdata[[binary_col]] <- ifelse(pdata[[binary_col]] == "High", 1, 0)
    }
    variable <- paste0(variable, "_binary")
  }

  # Run Cox models
  result_list <- lapply(
    variable,
    function(var) {
      fit <- survival::coxph(
        survival::Surv(pdata$time_iobr, pdata$status_iobr) ~ pdata[[var]],
        data = pdata
      )
      getHRandCIfromCoxph(fit)
    }
  )

  # Combine results
  result <- dplyr::bind_rows(result_list)
  result$ID <- variable

  if (best_cutoff) {
    result$ID <- gsub("_binary$", "", result$ID)
  }

  result <- result[order(result$P, decreasing = FALSE), ]
  rownames(result) <- NULL

  tibble::as_tibble(result)
}
