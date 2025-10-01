#' Extract best cut-off and add new binary object to data frame
#'
#' @description The "best_cutoff" function is used to find the best cutoff point for a continuous variable in survival analysis. It takes input data ("pdata") containing a continuous variable ("variable") and survival information ("time" and "status"). It returns the modified input data with a new binary variable created based on the best cutoff point.
#' @param pdata phenotype data with survival information and features
#' @param variable The name of the continuous variable in the input data for which the best cutoff needs to be determined.
#' @param time The name of the column in the input data representing the time-to-event (survival time). The default value is "time".
#' @param status The name of the column in the input data representing the event status (censoring information). The default value is "status".
#' @param PrintResult A logical value indicating whether to print the results. The default value is TRUE.
#'
#' @return pdata with binary variables and best cutoff
#' @export
#' @author Dongqiang Zeng
#' @examples
#'
#' # Loading TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Finding the best cutoff value of TMEscore for survival
#' sig_stad2 <- best_cutoff2(pdata = sig_stad, variable = "TMEscore_CIR", time = "OS_time", status = "OS_status", PrintResult = T)
#' sig_stad2$best_cutoff
#'
best_cutoff2 <- function(pdata, variable, time = "time", status = "status", PrintResult = T) {
  pdata <- as.data.frame(pdata)
  colnames(pdata)[which(colnames(pdata) == time)] <- "time"
  colnames(pdata)[which(colnames(pdata) == status)] <- "status"

  pdata <- pdata[!is.na(pdata$time), ]
  pdata <- pdata[!is.na(pdata$status), ]

  pdata$time <- as.numeric(pdata$time)
  pdata$status <- as.numeric(pdata$status)

  y <- Surv(pdata$time, pdata$status)
  iscutoff <- surv_cutpoint(pdata, time = "time", event = "status", variables = variable)
  aa <- paste("best cutoff = ", iscutoff$cutpoint$cutpoint)
  plot(iscutoff, variable, palette = "npg")

  bb <- base::summary(coxph(y ~ pdata[, which(colnames(pdata) == variable)], data = pdata))

  variable2 <- paste(variable, "_binary", sep = "")
  pdata[, variable2] <- ifelse(pdata[, variable] < iscutoff$cutpoint$cutpoint, "Low", "High")
  pdata[, variable2] <- as.factor(pdata[, variable2])

  colnames(pdata)[which(colnames(pdata) == "time")] <- time
  colnames(pdata)[which(colnames(pdata) == "status")] <- status

  cc <- summary(pdata[, variable2])
  dd <- base::summary(coxph(y ~ pdata[, which(colnames(pdata) == variable2)], data = pdata))
  if (PrintResult) {
    print(list(best_cutoff = aa, cox_continuous_object = bb, summary_binary_variable = cc, cox_binary_object = dd))
  }
  return(list(pdata = pdata, best_cutoff = iscutoff$cutpoint$cutpoint, cox_continuous_object = bb, summary_binary_variable = cc, cox_binary_object = dd))
}
