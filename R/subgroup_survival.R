#' Subgroup Survival Analysis Using Cox Proportional Hazards Models
#'
#' Extracts hazard ratios (HR) and 95% confidence intervals from Cox proportional hazards models across specified subgroups.
#'
#' @param pdata Data frame containing variables, follow-up time, and outcome.
#' @param variables Character vector. Subgrouping variables (each processed independently).
#' @param time_name Character. Column name of follow-up time. Default is "time".
#' @param status_name Character. Column name of event status (1/0). Default is "status".
#' @param object Character. Variable of interest used in Cox model. If it has levels "High" and "Low", recode "High" to 0 and "Low" to 1 before calling.
#'
#' @return Data frame summarizing subgroup Cox results (HR, CI, p-value).
#' @author Dongqiang Zeng
#' @export
#' @import survival
#' @examples
#' data(subgroup_data)
#' input <- subset(subgroup_data, time > 0 & !is.na(status) & !is.na(AJCC_stage))
#' # Binary variable example
#' res_bin <- subgroup_survival(pdata = input, time_name = "time", status_name = "status",
#'   variables = c("ProjectID", "AJCC_stage"), object = "score_binary")
#' # Continuous variable example
#' res_cont <- subgroup_survival(pdata = input, time_name = "time", status_name = "status",
#'   variables = c("ProjectID", "AJCC_stage"), object = "score")

subgroup_survival <- function(pdata, time_name = "time", status_name = "status", variables, object) {
  P <- 1
  HR <- 1
  CI_low_0.95 <- 1
  CI_up_0.95 <- 1
  result <- data.frame(P, HR, CI_low_0.95, CI_up_0.95, row.names = "defult")
  pdata <- as.data.frame(pdata)
  for (sig in variables) {
    ind <- which(colnames(pdata) == sig)

    tmp <- pdata[!is.na(pdata[, ind]), ]
    if (dim(tmp)[1] == 0) {
      next
    }
    tmp[, ind] <- as.factor(as.character(tmp[, ind]))

    result2 <- data.frame(P, HR, CI_low_0.95, CI_up_0.95, row.names = "defult")
    nl <- nlevels(tmp[, ind])
    for (i in 1:nl) {
      tmp1 <- tmp[as.character(tmp[, ind]) == names(summary(tmp[, ind]))[i], ]
      if (dim(tmp1)[1] == 0) {
        next
      }
      fit <- coxph(Surv(tmp1[, time_name], tmp1[, status_name]) ~ tmp1[, object], data = tmp1)
      result1 <- getHRandCIfromCoxph(fit)
      rownames(result1) <- paste(sig, levels(tmp[, ind])[i], sep = "_")
      result2 <- rbind(result2, result1)
    }
    result <- rbind(result, result2)
  }
  result[result > 100] <- Inf
  return(result[-c(grep(rownames(result), pattern = "defult")), ])
}
########################################
