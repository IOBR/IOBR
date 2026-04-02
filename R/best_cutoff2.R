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
#' set.seed(123)
#' pdata <- data.frame(
#'   time = rexp(100),
#'   status = rbinom(100, 1, 0.5),
#'   score = rnorm(100, mean = 50, sd = 10)
#' )
#' result <- best_cutoff2(pdata, variable = "score", print_result = FALSE)
#' result$best_cutoff
best_cutoff2 <- function(pdata, variable, time = "time", status = "status",
                         print_result = TRUE) {
  if (!is.data.frame(pdata)) {
    cli::cli_abort("{.arg pdata} must be a data frame")
  }
  if (nrow(pdata) == 0) {
    cli::cli_abort("{.arg pdata} has no rows")
  }
  if (!variable %in% colnames(pdata)) {
    cli::cli_abort("Variable {.val {variable}} not found in pdata")
  }
  if (!time %in% colnames(pdata)) {
    cli::cli_abort("Time column {.val {time}} not found in pdata")
  }
  if (!status %in% colnames(pdata)) {
    cli::cli_abort("Status column {.val {status}} not found in pdata")
  }

  pdata <- as.data.frame(pdata)
  colnames(pdata)[colnames(pdata) == time] <- "time"
  colnames(pdata)[colnames(pdata) == status] <- "status"

  pdata <- pdata[!is.na(pdata$time) & !is.na(pdata$status), , drop = FALSE]
  if (nrow(pdata) == 0) {
    cli::cli_abort("No valid data after removing NA values")
  }

  pdata$time <- as.numeric(pdata$time)
  pdata$status <- as.numeric(pdata$status)

  surv_obj <- survival::Surv(pdata$time, pdata$status)
  iscutoff <- survminer::surv_cutpoint(
    pdata,
    time = "time",
    event = "status",
    variables = variable
  )

  cutoff_value <- iscutoff$cutpoint$cutpoint
  cli::cli_alert_success("Best cutoff for {.val {variable}}: {round(cutoff_value, 3)}")
  plot(iscutoff, variable, palette = "npg")

  cox_cont <- summary(survival::coxph(
    surv_obj ~ pdata[[variable]],
    data = pdata
  ))

  variable2 <- paste0(variable, "_binary")
  pdata[[variable2]] <- ifelse(
    pdata[[variable]] < cutoff_value,
    "Low",
    "High"
  )
  pdata[[variable2]] <- factor(pdata[[variable2]], levels = c("Low", "High"))

  colnames(pdata)[colnames(pdata) == "time"] <- time
  colnames(pdata)[colnames(pdata) == "status"] <- status

  binary_summary <- summary(pdata[[variable2]])
  cox_binary <- summary(survival::coxph(
    surv_obj ~ pdata[[variable2]],
    data = pdata
  ))

  if (print_result) {
    print(list(
      best_cutoff = paste0("best cutoff = ", cutoff_value),
      cox_continuous_object = cox_cont,
      summary_binary_variable = binary_summary,
      cox_binary_object = cox_binary
    ))
  }

  list(
    pdata = pdata,
    best_cutoff = cutoff_value,
    cox_continuous_object = cox_cont,
    summary_binary_variable = binary_summary,
    cox_binary_object = cox_binary
  )
}
