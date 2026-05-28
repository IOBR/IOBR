#' Subgroup Survival Analysis Using Cox Proportional Hazards Models
#'
#' @description
#' Extracts hazard ratios (HR) and 95% confidence intervals from Cox
#' proportional hazards models across specified subgroups.
#'
#' @param pdata Data frame containing variables, follow-up time, and outcome.
#' @param variables Character vector. Subgrouping variables (each processed
#'   independently).
#' @param time_name Character. Column name of follow-up time. Default is `"time"`.
#' @param status_name Character. Column name of event status (1/0).
#'   Default is `"status"`.
#' @param object Character. Variable of interest used in Cox model. If it has
#'   levels "High" and "Low", recode "High" to 0 and "Low" to 1 before calling.
#'
#' @return Data frame summarizing subgroup Cox results (HR, CI, p-value).
#'
#' @export
#' @author Dongqiang Zeng
#' @import survival
#'
#' @examples
#' \donttest{
#' subgroup_data <- load_data("subgroup_data")
#' input <- subset(subgroup_data, time > 0 & !is.na(status) & !is.na(AJCC_stage))
#'
#' # Binary variable example
#' res_bin <- subgroup_survival(
#'   pdata = input, time_name = "time", status_name = "status",
#'   variables = c("ProjectID", "AJCC_stage"), object = "score_binary"
#' )
#'
#' # Continuous variable example
#' res_cont <- subgroup_survival(
#'   pdata = input, time_name = "time", status_name = "status",
#'   variables = c("ProjectID", "AJCC_stage"), object = "score"
#' )
#' }
subgroup_survival <- function(pdata,
                              time_name = "time",
                              status_name = "status",
                              variables,
                              object) {
  if (!is.data.frame(pdata)) {
    cli::cli_abort("{.arg pdata} must be a data frame.")
  }

  required_cols <- c(time_name, status_name, object, variables)
  missing_cols <- setdiff(required_cols, colnames(pdata))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing required column(s): {.val {missing_cols}}")
  }

  pdata <- as.data.frame(pdata)

  result <- data.frame(
    P = 1, HR = 1, CI_low_0.95 = 1, CI_up_0.95 = 1,
    row.names = "default"
  )

  for (sig in variables) {
    tmp <- pdata[!is.na(pdata[[sig]]), ]

    if (nrow(tmp) == 0) {
      cli::cli_warn("No valid data for variable {.val {sig}}, skipping.")
      next
    }

    tmp[[sig]] <- as.factor(as.character(tmp[[sig]]))

    result2 <- data.frame(
      P = 1, HR = 1, CI_low_0.95 = 1, CI_up_0.95 = 1,
      row.names = "default"
    )

    nl <- nlevels(tmp[[sig]])

    for (i in seq_len(nl)) {
      level_name <- levels(tmp[[sig]])[i]
      tmp1 <- tmp[as.character(tmp[[sig]]) == level_name, ]

      if (nrow(tmp1) == 0) {
        next
      }

      fit <- survival::coxph(
        survival::Surv(tmp1[[time_name]], tmp1[[status_name]]) ~ tmp1[[object]],
        data = tmp1
      )
      result1 <- getHRandCIfromCoxph(fit)
      rownames(result1) <- paste(sig, level_name, sep = "_")
      result2 <- rbind(result2, result1)
    }

    result <- rbind(result, result2)
  }

  result[result > 100] <- Inf
  result <- result[!grepl("default", rownames(result)), , drop = FALSE]

  result
}
