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
#'   Default is \code{"time"}.
#' @param status Character string specifying the column name containing event status
#'   (1 = event occurred, 0 = censored). Default is \code{"status"}.
#' @param best_cutoff Logical indicating whether to compute optimal cutoffs for
#'   continuous variables and analyze dichotomized versions. Default is \code{FALSE}.
#'
#' @return Data frame containing hazard ratios (HR), 95% confidence intervals (CI),
#'   and p-values for each variable, sorted by p-value.
#'
#' @author Dongqiang Zeng
#' @export
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Perform batch survival analysis
#' batch_surv(pdata = sig_stad, variable = colnames(sig_stad)[69:ncol(sig_stad)],
#'            time = "OS_time", status = "OS_status")
batch_surv <- function(pdata, variable, time = "time", status = "status", best_cutoff = FALSE) {
  pdata <- as.data.frame(pdata)

  colnames(pdata)[which(colnames(pdata) == time)] <- "time"
  colnames(pdata)[which(colnames(pdata) == status)] <- "status"

  pdata <- pdata[!is.na(pdata$time), ]
  pdata <- pdata[!is.na(pdata$status), ]

  pdata$time <- as.numeric(pdata$time)
  pdata$status <- as.numeric(pdata$status)
  #################################################
  if (best_cutoff == TRUE) {
    for (i in 1:length(variable)) {
      var <- variable[i]
      pdata <- best_cutoff(pdata = pdata, time = "time", status = "status", variable = var, PrintResult = FALSE)
      pdata[, paste0(var, "_binary")] <- ifelse(pdata[, paste0(var, "_binary")] == "High", 1, 0)
    }
    # print(pdata)
    variable <- paste0(variable, "_binary")
  }
  #################################################
  result_list <- lapply(pdata[, variable], function(x) coxph(Surv(pdata$time, pdata$status) ~ x, data = pdata[, variable]))
  result <- data.frame(NULL)
  for (i in 1:length(result_list)) {
    result1 <- getHRandCIfromCoxph(result_list[[i]])
    rownames(result1) <- variable[i]
    result <- rbind(result, result1)
  }
  # result[result>1000]<-Inf

  result <- result[order(result$P, decreasing = F), ]
  result <- tibble::rownames_to_column(result, var = "ID")

  if (best_cutoff == TRUE) {
    result$ID <- gsub(result$ID, pattern = "_binary", replacement = "")
  }
  result <- tibble::as_tibble(result)
  return(result)
}
