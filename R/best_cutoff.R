#' Extract Best Cutoff and Add Binary Variable to Data Frame
#'
#' @description
#' Determines the optimal cutoff point for a continuous variable in survival analysis
#' using the maximally selected rank statistics method. Creates a binary variable based
#' on the identified cutoff and adds it to the input data frame.
#'
#' @param pdata Data frame containing survival information and the continuous variable.
#' @param variable Character string specifying the name of the continuous variable for
#'   which the optimal cutoff should be determined.
#' @param time Character string specifying the column name containing time-to-event data.
#'   Default is \code{"time"}.
#' @param status Character string specifying the column name containing event status
#'   (censoring information). Default is \code{"status"}.
#' @param PrintResult Logical indicating whether to print detailed results including
#'   cutoff value, Cox model summaries for continuous and binary variables. Default is
#'   \code{TRUE}.
#'
#' @return Data frame identical to \code{pdata} with an additional binary column named
#'   \code{<variable>_binary} containing "High" and "Low" categories based on the
#'   optimal cutoff.
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Find the best cutoff value of TMEscore for survival analysis
#' sig_stad2 <- best_cutoff(pdata = sig_stad, variable = "TMEscore_CIR",
#'                          time = "OS_time", status = "OS_status", PrintResult = TRUE)
#' table(sig_stad2$TMEscore_CIR_binary)
best_cutoff <- function(pdata, variable, time = "time", status = "status", PrintResult = T) {
  pdata <- as.data.frame(pdata)
  colnames(pdata)[which(colnames(pdata) == time)] <- "time"
  colnames(pdata)[which(colnames(pdata) == status)] <- "status"

  pdata <- pdata[!is.na(pdata$time), ]
  pdata <- pdata[!is.na(pdata$status), ]

  pdata$time <- as.numeric(pdata$time)
  pdata$status <- as.numeric(pdata$status)

  y <- Surv(pdata$time, pdata$status)
  iscutoff <- surv_cutpoint(pdata, time = "time", event = "status", variables = variable)

  aa <- paste(">>>-- The best cutoff is = ", iscutoff$cutpoint$cutpoint)
  message(aa)

  # plot(iscutoff,variable,palette="npg")

  bb <- base::summary(coxph(y ~ pdata[, which(colnames(pdata) == variable)], data = pdata))
  variable2 <- paste(variable, "_binary", sep = "")
  pdata[, variable2] <- ifelse(pdata[, variable] <= iscutoff$cutpoint$cutpoint, "Low", "High")
  pdata[, variable2] <- as.factor(pdata[, variable2])
  cc <- summary(pdata[, variable2])
  dd <- base::summary(coxph(y ~ pdata[, which(colnames(pdata) == variable2)], data = pdata))

  if (PrintResult) {
    print(list(best_cutoff = aa, cox_continuous_object = bb, summary_binary_variable = cc, cox_binary_object = dd))
  }

  colnames(pdata)[which(colnames(pdata) == "time")] <- time
  colnames(pdata)[which(colnames(pdata) == "status")] <- status
  return(pdata)
}
