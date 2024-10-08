##' Extract hazard ratio and confidence intervals from a coxph object
##'
##' This convenience function extracts hazard ratio and confidence intervals from a coxph object,
##' as generated by a call to coxph(), comparing survival times between two groups.
##'
##' @param coxphData coxph object; example usage:
##' \dontrun{
##' library(survival)
##' # Create example data
##' set.seed(123)
##' df <- data.frame(
##'   TTE = rexp(200, rate = 0.1),  # Time to event data
##'   Cens = rbinom(200, size = 1, prob = 0.7),  # Censoring information
##'   group = sample(c("Treatment", "Control"), 200, replace = TRUE)  # Group information
##' )
##' # Fit Cox proportional hazards model
##' coxphData <- coxph(Surv(TTE, Cens) ~ group, data = df)
##' # Extract HR and CI
##' results <- getHRandCIfromCoxph(coxphData)
##' print(results)
##' }
##'
##' @author Dorothee Nickles
##' @author Dongqiang Zeng
##' @export
##' @import survival
getHRandCIfromCoxph <- function(coxphData) {

  stopifnot(is(coxphData, "coxph"))

  tmp <- cbind(
    summary(coxphData)$coef[,
                            c("Pr(>|z|)", "exp(coef)"),
                            drop=FALSE],
    summary(coxphData)$conf[,
                            c("lower .95", "upper .95"),
                            drop=FALSE]
  )
  colnames(tmp) <- c("P","HR","CI_low_0.95","CI_up_0.95")
  tmp[,2:4]<-round(tmp[,2:4],digits=4)
  return(tmp)
}
