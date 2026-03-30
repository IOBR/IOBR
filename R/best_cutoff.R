#' Extract Best Cutoff and Add Binary Variable to Data Frame
#'
#' @description
#' Determines the optimal cutoff point for a continuous variable in survival
#' analysis using the maximally selected rank statistics method. Creates a
#' binary variable based on the identified cutoff and adds it to the input data
#' frame.
#'
#' @param pdata Data frame containing survival information and the continuous
#'   variable.
#' @param variable Character string specifying the name of the continuous
#'   variable for which the optimal cutoff should be determined.
#' @param time Character string specifying the column name containing
#'   time-to-event data. Default is `"time"`.
#' @param status Character string specifying the column name containing event
#'   status (censoring information). Default is `"status"`.
#' @param PrintResult Logical indicating whether to print detailed results
#'   including cutoff value and Cox model summaries. Default is `TRUE`.
#'
#' @return Data frame identical to `pdata` with an additional binary column
#'   named `<variable>_binary` containing "High" and "Low" categories based on
#'   the optimal cutoff.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' sig_stad <- load_data("sig_stad")
#' sig_stad2 <- best_cutoff(
#'   pdata = sig_stad,
#'   variable = "TMEscore_CIR",
#'   time = "OS_time",
#'   status = "OS_status",
#'   PrintResult = TRUE
#' )
#' table(sig_stad2$TMEscore_CIR_binary)
#' }
best_cutoff <- function(pdata, variable, time = "time",
                        status = "status", print_result = TRUE) {
  # Input validation
  if (!is.data.frame(pdata)) {
    cli::cli_abort("{.arg pdata} must be a data frame")
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

  # Suppress NOTES about unused variable
  surv_obj <- NULL

  pdata <- as.data.frame(pdata)

  # Rename columns temporarily
  colnames(pdata)[colnames(pdata) == time] <- "time_iobr"
  colnames(pdata)[colnames(pdata) == status] <- "status_iobr"

  # Remove NA values
  pdata <- pdata[!is.na(pdata$time_iobr) & !is.na(pdata$status_iobr), ]

  if (nrow(pdata) == 0) {
    cli::cli_abort("No valid data after removing NA values")
  }

  # Ensure numeric
  pdata$time_iobr <- as.numeric(pdata$time_iobr)
  pdata$status_iobr <- as.numeric(pdata$status_iobr)

  # Find best cutoff
  surv_obj <- survival::Surv(pdata$time_iobr, pdata$status_iobr)

  iscutoff <- survminer::surv_cutpoint(
    pdata,
    time = "time_iobr",
    event = "status_iobr",
    variables = variable
  )

  cutoff_value <- iscutoff$cutpoint$cutpoint
  cli::cli_alert_success(paste(
    "Best cutoff for {.val {variable}}:",
    "{round(cutoff_value, 3)}"
  ))

  # Cox model for continuous variable
  cox_cont <- summary(survival::coxph(
    surv_obj ~ pdata[[variable]],
    data = pdata
  ))

  # Create binary variable
  variable2 <- paste0(variable, "_binary")
  pdata[[variable2]] <- ifelse(
    pdata[[variable]] <= cutoff_value,
    "Low",
    "High"
  )
  pdata[[variable2]] <- factor(pdata[[variable2]], levels = c("Low", "High"))

  # Summary statistics
  binary_summary <- summary(pdata[[variable2]])

  # Cox model for binary variable
  cox_binary <- summary(survival::coxph(
    surv_obj ~ pdata[[variable2]],
    data = pdata
  ))

  if (print_result) {
    print(list(
      best_cutoff = paste0("The best cutoff is = ", round(cutoff_value, 3)),
      cox_continuous_object = cox_cont,
      summary_binary_variable = binary_summary,
      cox_binary_object = cox_binary
    ))
  }

  # Restore original column names
  colnames(pdata)[colnames(pdata) == "time_iobr"] <- time
  colnames(pdata)[colnames(pdata) == "status_iobr"] <- status

  pdata
}
