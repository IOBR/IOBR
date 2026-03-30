#' Extract Best Cutoff and Add Binary Variable to Data Frame
#'
#' @description
#' Finds the best cutoff point for a continuous variable in survival analysis.
#' Takes input data containing a continuous variable and survival information
#' (time and status). Returns the modified input data with a new binary variable
#' created based on the best cutoff point.
#'
#' @param pdata Data frame containing survival information and features.
#' @param variable Character string specifying the continuous variable name.
#' @param time Character string specifying the time-to-event column name.
#'   Default is `"time"`.
#' @param status Character string specifying the event status column name.
#'   Default is `"status"`.
#' @param print_result Logical indicating whether to print results.
#'   Default is `TRUE`.
#'
#' @return List containing:
#' \describe{
#'   \item{pdata}{Data frame with binary variable added}
#'   \item{best_cutoff}{Numeric cutoff value}
#'   \item{cox_continuous_object}{Cox model summary for continuous variable}
#'   \item{summary_binary_variable}{Summary of binary variable}
#'   \item{cox_binary_object}{Cox model summary for binary variable}
#' }
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' # Loading TCGA-STAD microenvironment signature data
#' sig_stad <- load_data("sig_stad")
#' # Finding the best cutoff value of TMEscore for survival
#' result <- best_cutoff2(
#'   pdata = sig_stad,
#'   variable = "TMEscore_CIR",
#'   time = "OS_time",
#'   status = "OS_status",
#'   print_result = TRUE
#' )
#' result$best_cutoff
#' }
best_cutoff2 <- function(pdata, variable, time = "time", status = "status",
                         print_result = TRUE) {
  pdata <- as.data.frame(pdata)
  colnames(pdata)[which(colnames(pdata) == time)] <- "time"
  colnames(pdata)[which(colnames(pdata) == status)] <- "status"

  pdata <- pdata[!is.na(pdata$time), ]
  pdata <- pdata[!is.na(pdata$status), ]

  pdata$time <- as.numeric(pdata$time)
  pdata$status <- as.numeric(pdata$status)

  surv_obj <- survival::Surv(pdata$time, pdata$status)
  iscutoff <- survminer::surv_cutpoint(
    pdata,
    time = "time",
    event = "status",
    variables = variable
  )
  aa <- paste("best cutoff = ", iscutoff$cutpoint$cutpoint)
  plot(iscutoff, variable, palette = "npg")

  bb <- base::summary(survival::coxph(
    surv_obj ~ pdata[, which(colnames(pdata) == variable)],
    data = pdata
  ))

  variable2 <- paste(variable, "_binary", sep = "")
  pdata[, variable2] <- ifelse(
    pdata[, variable] < iscutoff$cutpoint$cutpoint,
    "Low",
    "High"
  )
  pdata[, variable2] <- as.factor(pdata[, variable2])

  colnames(pdata)[which(colnames(pdata) == "time")] <- time
  colnames(pdata)[which(colnames(pdata) == "status")] <- status

  cc <- summary(pdata[, variable2])
  dd <- base::summary(survival::coxph(
    surv_obj ~ pdata[, which(colnames(pdata) == variable2)],
    data = pdata
  ))
  if (print_result) {
    print(list(
      best_cutoff = aa,
      cox_continuous_object = bb,
      summary_binary_variable = cc,
      cox_binary_object = dd
    ))
  }
  list(
    pdata = pdata,
    best_cutoff = iscutoff$cutpoint$cutpoint,
    cox_continuous_object = bb,
    summary_binary_variable = cc,
    cox_binary_object = dd
  )
}
