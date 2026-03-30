#' Add Risk Score to Dataset
#'
#' @description
#' Computes a risk score for each observation based on Cox proportional hazards
#' regression or binary logistic regression. The function fits the specified model
#' and returns the dataset with an added risk score column.
#'
#' @param input Data frame containing the variables for analysis.
#' @param family Character string specifying the model family: `"cox"` for Cox
#'   proportional hazards regression or `"binary"` for logistic regression.
#'   Default is `"cox"`.
#' @param target Character string specifying the target variable name. Required when
#'   `family = "binary"`.
#' @param time Character string specifying the time-to-event variable name. Required
#'   when `family = "cox"`.
#' @param status Character string specifying the event status variable name. Required
#'   when `family = "cox"`.
#' @param vars Character vector of variable names to include in the model.
#' @param new_var_name Character string specifying the name for the new risk score
#'   column. Default is `"riskscore"`.
#'
#' @return Data frame identical to `input` with an additional column containing
#'   risk scores (linear predictors for Cox models or predicted probabilities for
#'   logistic models).
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   result <- add_riskscore(
#'     lung,
#'     time = "time", status = "status",
#'     vars = c("age", "sex")
#'   )
#'   head(result)
#' }
#' }
add_riskscore <- function(input, family = c("cox", "binary"), target = NULL,
                          time = NULL, status = NULL, vars,
                          new_var_name = "riskscore") {

  family <- rlang::arg_match(family)

  # Input validation
  if (!is.data.frame(input)) {
    cli::cli_abort("{.arg input} must be a data frame")
  }
  if (nrow(input) == 0) {
    cli::cli_abort("{.arg input} has no rows")
  }
  if (!is.character(vars) || length(vars) == 0) {
    cli::cli_abort("{.arg vars} must be a non-empty character vector")
  }

  missing_vars <- setdiff(vars, colnames(input))
  if (length(missing_vars) > 0) {
    cli::cli_abort("Variables not found in input: {.val {missing_vars}}")
  }

  # Validate family-specific parameters
  if (family == "cox") {
    if (is.null(time) || !time %in% colnames(input)) {
      cli::cli_abort("Cox model requires {.arg time} to be a valid column name")
    }
    if (is.null(status) || !status %in% colnames(input)) {
      cli::cli_abort("Cox model requires {.arg status} to be a valid column name")
    }
  } else if (family == "binary") {
    if (is.null(target) || !target %in% colnames(input)) {
      cli::cli_abort("Binary model requires {.arg target} to be a valid column name")
    }
    if (length(unique(input[[target]])) != 2) {
      cli::cli_abort("Target variable must have exactly 2 unique values")
    }
  }

  # Build and fit model
  if (family == "cox") {
    formula <- stats::as.formula(
      paste("survival::Surv(", time, ",", status, ") ~", paste(vars, collapse = " + "))
    )
    model <- survival::coxph(formula, data = input)
    input[[new_var_name]] <- stats::predict(model, newdata = input, type = "lp")
  } else {
    input[[target]] <- as.numeric(factor(input[[target]])) - 1
    formula <- stats::as.formula(
      paste(target, "~", paste(vars, collapse = " + "))
    )
    model <- stats::glm(formula, data = input, family = stats::binomial())
    input[[new_var_name]] <- stats::predict(model, newdata = input, type = "response")
  }

  if (interactive()) {
    print(summary(model))
  }

  input
}
