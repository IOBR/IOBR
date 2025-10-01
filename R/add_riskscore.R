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
#' # Example with survival data
#' data(ovarian)
#' ovarian$rscore <- add_riskscore(ovarian, time = "time", status = "status",
#'                                  vars = c("resid.ds", "rx", "ecog.ps"))
add_riskscore <- function(input, family = "cox", target = NULL, time = NULL, status = NULL, vars, new_var_name = "riskscore") {
  # Check the 'family' parameter and perform the corresponding model fitting
  if (family == "cox") {
    # Combine all variables to be used in the survival model
    formula <- as.formula(paste("Surv(", time, ",", status, ") ~ ", paste(vars, collapse = " + ")))

    # Fit the Cox survival model
    model <- coxph(formula, data = input)

    # Print model summary
    print(summary(model))

    # Add a risk score column to the input data frame
    input[[new_var_name]] <- predict(model, newdata = input, type = "lp")
  } else if (family == "binary") {
    # Ensure the target is binary
    if (length(unique(input[[target]])) != 2) {
      stop("Target variable must be binary.")
    }

    # Convert the target variable to 0 and 1
    input[[target]] <- as.numeric(as.factor(input[[target]])) - 1

    # Combine all variables to be used in the logistic regression
    formula <- as.formula(paste(target, "~", paste(vars, collapse = " + ")))

    # Fit the logistic regression model
    model <- glm(formula, data = input, family = binomial())

    # Print model summary
    print(summary(model))

    # Add a prediction column to the input data frame
    input[[new_var_name]] <- predict(model, newdata = input, type = "response")
  } else {
    stop("Unsupported family specified. Use 'cox' for Cox regression or 'binary' for binary logistic regression.")
  }

  # Return the input data frame containing the risk scores or predictions
  return(input)
}
