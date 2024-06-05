#' Add Risk Score to Dataset
#'
#' This function computes a risk score for each observation in the dataset
#' based on the provided Cox proportional hazards model or logistic regression.
#' @param input A data frame containing the data to be analyzed.
#' @param family A character string specifying whether to use "cox" for Cox regression
#'        or "binary" for logistic regression. Defaults to "cox".
#' @param target The name of the target variable (only used if family = "binary").
#' @param time The name of the time to event variable (only used if family = "cox").
#' @param status The name of the event status variable (only used if family = "cox").
#' @param vars A vector of variable names to include in the model.
#' @param new_var_name The name of the new variable to be created for storing risk scores.
#'        Defaults to "riskscore".
#' @return The input data frame with an additional column containing the risk scores.
#' @author Dongqiang Zeng
#' @examples
#' data(ovarian)
#' ovarian$rscore <- add_riskscore(ovarian, time = "time", status = "status", vars=c("resid.ds", "rx", "ecog.ps"))

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
