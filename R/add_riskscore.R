#' Add Risk Score to Dataset
#'
#' @description
#' Computes a risk score for each observation based on Cox proportional hazards
#' regression or binary logistic regression. The function fits the specified model
#' and returns the dataset with an added risk score column.
#'
#' @param input Data frame containing the variables for analysis.
#' @param family Character string specifying the model family: \code{"cox"} for Cox
#'   proportional hazards regression or \code{"binary"} for logistic regression.
#'   Default is \code{"cox"}.
#' @param target Character string specifying the target variable name. Required when
#'   \code{family = "binary"}.
#' @param time Character string specifying the time-to-event variable name. Required
#'   when \code{family = "cox"}.
#' @param status Character string specifying the event status variable name. Required
#'   when \code{family = "cox"}.
#' @param vars Character vector of variable names to include in the model.
#' @param new_var_name Character string specifying the name for the new risk score
#'   column. Default is \code{"riskscore"}.
#'
#' @return Data frame identical to \code{input} with an additional column containing
#'   risk scores (linear predictors for Cox models or predicted probabilities for
#'   logistic models).
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$rscore <- add_riskscore(
#'     lung,
#'     time = "time", status = "status",
#'     vars = c("age", "sex")
#'   )
#' }
add_riskscore <- function(input, family = "cox", target = NULL, time = NULL, status = NULL, vars, new_var_name = "riskscore") {
  # Input validation
  if (is.null(input) || !is.data.frame(input)) {
    stop("'input' must be a non-null data frame.")
  }
  if (nrow(input) == 0) {
    stop("'input' has no rows.")
  }
  if (!family %in% c("cox", "binary")) {
    stop("'family' must be either 'cox' or 'binary'.")
  }
  if (missing(vars) || !is.character(vars) || length(vars) == 0) {
    stop("'vars' must be a non-empty character vector of variable names.")
  }
  missing_vars <- setdiff(vars, colnames(input))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variables not found in input: %s", paste(missing_vars, collapse = ", ")))
  }

  # Validate family-specific parameters
  if (family == "cox") {
    if (is.null(time) || !time %in% colnames(input)) {
      stop(sprintf("Cox model requires 'time' parameter pointing to a valid column."))
    }
    if (is.null(status) || !status %in% colnames(input)) {
      stop(sprintf("Cox model requires 'status' parameter pointing to a valid column."))
    }
  } else if (family == "binary") {
    if (is.null(target) || !target %in% colnames(input)) {
      stop(sprintf("Binary model requires 'target' parameter pointing to a valid column."))
    }
    if (length(unique(input[[target]])) != 2) {
      stop("Target variable must be binary (exactly 2 unique values).")
    }
  }

  # Build and fit model
  if (family == "cox") {
    formula <- stats::as.formula(paste("survival::Surv(", time, ",", status, ") ~ ", paste(vars, collapse = " + ")))
    model <- survival::coxph(formula, data = input)
    input[[new_var_name]] <- predict(model, newdata = input, type = "lp")
  } else {
    # Convert target to 0/1 for binary logistic
    input[[target]] <- as.numeric(factor(input[[target]])) - 1
    formula <- stats::as.formula(paste(target, "~", paste(vars, collapse = " + ")))
    model <- stats::glm(formula, data = input, family = stats::binomial())
    input[[new_var_name]] <- predict(model, newdata = input, type = "response")
  }

  if (interactive()) {
    print(summary(model))
  }

  return(input)
}
