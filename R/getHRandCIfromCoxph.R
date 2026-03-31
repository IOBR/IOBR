#' Extract Hazard Ratio and Confidence Intervals from Cox Model
#'
#' @description
#' Extracts hazard ratio (HR) and 95% confidence intervals from a fitted
#' Cox proportional hazards model.
#'
#' @param coxphData Fitted Cox model object from `survival::coxph()`.
#'
#' @return Data frame with p-values, HR, and confidence intervals.
#'
#' @export
#' @author Dorothee Nickles
#' @author Dongqiang Zeng
#'
#' @examples
#' library(survival)
#' set.seed(123)
#' df <- data.frame(
#'   TTE = rexp(200, rate = 0.1),
#'   Cens = rbinom(200, size = 1, prob = 0.7),
#'   group = sample(c("Treatment", "Control"), 200, replace = TRUE)
#' )
#' coxphData <- survival::coxph(survival::Surv(TTE, Cens) ~ group, data = df)
#' results <- getHRandCIfromCoxph(coxphData)
#' print(results)
getHRandCIfromCoxph <- function(coxphData) {
  if (!inherits(coxphData, "coxph")) {
    cli::cli_abort("{.arg coxphData} must be a coxph object")
  }

  tmp <- cbind(
    summary(coxphData)$coef[, c("Pr(>|z|)", "exp(coef)"), drop = FALSE],
    summary(coxphData)$conf[, c("lower .95", "upper .95"), drop = FALSE]
  )
  colnames(tmp) <- c("P", "HR", "CI_low_0.95", "CI_up_0.95")
  tmp[, 2:4] <- round(tmp[, 2:4], digits = 4)

  as.data.frame(tmp)
}
