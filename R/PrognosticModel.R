#' Build Prognostic Models Using LASSO and Ridge Regression
#'
#' This function prepares data, splits it into training and testing sets, and fits LASSO and Ridge regression models
#' for survival analysis. It evaluates the models and optionally plots time-dependent ROC curves.
#'
#' @param x A matrix or data frame of predictor variables.
#' @param y A data frame of survival outcomes, including time and status columns.
#' @param scale Logical indicating whether to scale predictor variables. Default is FALSE.
#' @param seed Integer seed for random number generation to ensure reproducibility.
#' @param train_ratio Numeric proportion of data for training (e.g., 0.7). Default is 0.7.
#' @param nfold Integer number of folds for cross-validation. Default is 10.
#' @param plot Logical indicating whether to plot ROC curves. Default is TRUE.
#' @param cols Optional vector of colors for ROC curves. If NULL, uses default palette.
#' @param palette String specifying color palette. Default is "jama".
#'
#' @return A list containing results from LASSO and Ridge models, training data, and optionally ROC plots.
#' @importFrom glmnet cv.glmnet
#' @import dplyr
#' @import ggplot2
#' @export
#' @examples
#' data("imvigor210_sig", package = "IOBR")
#' data("imvigor210_pdata", package = "IOBR")
#' pdata_prog <- imvigor210_pdata %>%
#'   dplyr::select(ID, OS_days, OS_status) %>%
#'   mutate(OS_days = as.numeric(OS_days), OS_status = as.numeric(OS_status))
#' prognostic_result <- PrognosticModel(x = imvigor210_sig, y = pdata_prog,
#'                                      scale = TRUE, seed = 123456,
#'                                      train_ratio = 0.7, nfold = 10, plot = TRUE)
PrognosticModel <- function(x, y, scale = FALSE, seed = 123456, train_ratio = 0.7, nfold = 10, plot = TRUE, palette = "jama", cols = NULL) {
  x <- as.data.frame(x)
  y <- as.data.frame(y)

  message(paste0(">>> Processing data"))
  processdat <- ProcessingData(x = x, y = y, scale = scale, type = "survival")
  x_scale <- processdat$x_scale
  y <- processdat$y
  x_ID <- processdat$x_ID

  message(paste0(">>> Spliting data into training and validation data"))
  train_test <- SplitTrainTest(
    x = x_scale, y = y, train_ratio = train_ratio,
    type = "survival", seed = seed
  )
  train.x <- train_test$train.x
  train.y <- train_test$train.y
  test.x <- train_test$test.x
  test.y <- train_test$test.y
  train_sample <- train_test$train_sample
  return.x <- data.frame(ID = x_ID[train_sample], train.x)

  message(paste0(">>> Running ", "LASSO"))
  set.seed(seed)
  lasso_model <- glmnet::cv.glmnet(
    x = train.x, y = as.matrix(train.y),
    family = "cox", alpha = 1, nfolds = nfold
  )
  lasso_result <- PrognosticResult(model = lasso_model, train.x, train.y, test.x, test.y)
  if (plot) {
    p1 <- PlotTimeROC(
      train.x = train.x, train.y = train.y,
      test.x = test.x, test.y = test.y, model = lasso_model, modelname = "LASSO"
    )
  }

  message(paste0(">>> Running ", "RIDGE REGRESSION"))

  set.seed(seed)
  ridge_model <- glmnet::cv.glmnet(x = train.x, y = as.matrix(train.y), family = "cox", alpha = 0, nfolds = nfold)
  ridge_result <- PrognosticResult(model = ridge_model, train.x, train.y, test.x, test.y)

  if (plot) {
    p2 <- PlotTimeROC(
      train.x = train.x, train.y = train.y, cols = cols, palette = palette,
      test.x = test.x, test.y = test.y, model = ridge_model, modelname = "RIDGE"
    )
  }
  message(paste0(">>> Done !"))
  p <- p1 + p2
  print(p)
  return(list(lasso_result = lasso_result, ridge_result = ridge_result, train.x = return.x))
}

#' Prognostic Results for Survival Models
#'
#' This function computes and compiles the prognostic results from a survival model, specifically designed for models
#' like those created with `glmnet`. It calculates model coefficients at specific lambda values (`lambda.min` and `lambda.1se`)
#' and computes Area Under the Curve (AUC) metrics for both training and testing datasets. The AUC calculations
#' are performed at two regularization strengths for both train and test datasets.
#'
#' @param model A survival model object from which predictions and coefficients will be extracted.
#' @param train.x Matrix or data frame of training predictors that were used to fit the model.
#' @param train.y Training dataset outcomes, including survival time and status.
#' @param test.x Matrix or data frame of testing predictors for evaluating the model.
#' @param test.y Testing dataset outcomes, including survival time and status.
#'
#' @return A list containing three elements:
#'   - `model`: The model object passed to the function.
#'   - `coefs`: A data frame of model coefficients extracted for `lambda.min` and `lambda.1se`.
#'   - `AUC`: A data frame with the AUC values for training and testing datasets at both lambda values.
#'
#' @importFrom stats coef
#' @importFrom purrr pmap
#' @examples
#' # Assuming 'fit' is a Cox model fitted using `glmnet`
#' train_data <- list(x = matrix(rnorm(100 * 10), ncol = 10), y = Surv(rexp(100), rbinom(100, 1, 0.5)))
#' test_data <- list(x = matrix(rnorm(100 * 10), ncol = 10), y = Surv(rexp(100), rbinom(100, 1, 0.5)))
#' results <- PrognosticResult(
#'   model = fit, train.x = train_data$x, train.y = train_data$y,
#'   test.x = test_data$x, test.y = test_data$y
#' )
#' @export
PrognosticResult <- function(model, train.x, train.y, test.x, test.y) {
  # Extract coefficients at specified lambda values
  coefs <- cbind(stats::coef(model, s = "lambda.min"), stats::coef(model, s = "lambda.1se"))
  coefs <- data.frame(feature = rownames(coefs), lambda.min = coefs[, 1], lambda.1se = coefs[, 2])

  # Prepare data for AUC calculation
  newx <- list(train.x, train.x, test.x, test.x)
  s <- list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y <- list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)

  # Calculate AUC using the PrognosticAUC function
  auc <- args %>%
    purrr::pmap(PrognosticAUC, model = model) %>%
    do.call(rbind, .)
  rownames(auc) <- c("Train_lambda.min", "Train_lambda.1se", "Test_lambda.min", "Test_lambda.1se")

  # Return results as a list
  resultreturn <- list(model = model, coefs = coefs, AUC = auc)
}




#' Calculate Prognostic Area Under the Curve (AUC)
#'
#' This function evaluates the prognostic ability of a survival model by calculating the Area Under the
#' Curve (AUC) of time-dependent Receiver Operating Characteristic (ROC) curves at specified time points.
#' The function uses predictions made by the model to compute ROC statistics and AUC at the 30th and 90th
#' percentiles of survival time, which are commonly used thresholds to assess short-term and long-term risk.
#'
#' @param model A survival model object from which predictions will be made.
#'        This model should be capable of generating risk scores, such as a Cox proportional hazards model.
#' @param newx A matrix or data frame containing new input data for which predictions are to be made.
#'        This data should have the same features as used to train the model.
#' @param s The value of the penalty parameter 'lambda' at which predictions are requested.
#'        This can be a specific value or a character string specifying a criterion,
#'        such as 'lambda.min' or 'lambda.1se', commonly used in models like those from `glmnet`.
#' @param acture.y A data frame containing actual survival data, expected to have at least two columns:
#'        'time', which contains the survival time, and 'status', which is the event indicator (1 if the event occurred, 0 otherwise).
#'
#' @return A data frame with two columns containing the AUC values for the 30th and 90th percentile survival times.
#'         These are named 'probs.3' and 'probs.9' respectively.
#'
#' @importFrom stats predict
#' @importFrom timeROC timeROC
#' @examples
#' # Assuming 'fit' is a fitted survival model such as from `coxph` or `glmnet`
#' new_data <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' actual_outcome <- data.frame(time = rexp(100, rate = 0.1), status = rbinom(100, size = 1, prob = 0.5))
#' auc_results <- PrognosticAUC(fit, newx = new_data, s = "lambda.min", acture.y = actual_outcome)
#' @export
PrognosticAUC <- function(model, newx, s, acture.y) {
  riskscore <- stats::predict(model, newx = newx, s = s)
  timerocDat <- data.frame(risk = riskscore[, 1], acture.y)
  ROC <- with(
    timerocDat,
    timeROC::timeROC(
      T = time, delta = status,
      marker = risk, cause = 1,
      weighting = "marginal",
      time = quantile(time, probs = c(0.3, 0.9)),
      ROC = TRUE,
      iid = TRUE
    )
  )
  AUC <- data.frame(probs.3 = ROC$AUC[1], probs.9 = ROC$AUC[2])
  return(AUC)
}



#' Calculate Time-Dependent ROC
#'
#' This function computes the time-dependent Receiver Operating Characteristic (ROC) curve for a given survival model.
#' It uses the `timeROC` package to compute ROC statistics based on predicted risk scores from the model and actual survival outcomes.
#' This is particularly useful for evaluating the predictive accuracy of survival models at a specified time probability.
#'
#' @param model A survival model object from which predictions will be made.
#' @param newx A matrix or data frame of new data points for which predictions need to be made.
#' @param s The value of the penalty parameter 'lambda' at which predictions are requested.
#'         This can be a specific value or a character string specifying a criterion,
#'         such as 'lambda.min' or 'lambda.1se' which are commonly used in models like those from `glmnet`.
#' @param acture.y A data frame containing the actual survival data, expected to have at least two columns:
#'         'time' which is the survival time, and 'status' which is the event indicator (1 if the event occurred, 0 otherwise).
#' @param modelname A character string specifying the name of the model, used for identification in outputs or plots.
#' @param time_prob A numeric value specifying the quantile of interest for the ROC calculation.
#'         The default is 0.9, which calculates the ROC at the 90th percentile of the survival time.
#'
#' @return An object of class `timeROC` that contains the time-dependent ROC curve information.
#' @importFrom stats predict
#' @importFrom timeROC timeROC
#' @examples
#' # Assuming 'fit' is a Cox proportional hazards model fitted using `coxph` or `glmnet`
#' new_data <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' actual_outcome <- data.frame(time = rexp(100, rate = 0.1), status = rbinom(100, size = 1, prob = 0.5))
#' roc_info <- CalculateTimeROC(fit, newx = new_data, s = "lambda.min", acture.y = actual_outcome, modelname = "Cox Model")
#' @export
CalculateTimeROC <- function(model, newx, s, acture.y, modelname, time_prob = 0.9) {
  riskscore <- stats::predict(model, newx = newx, s = s)
  timerocDat <- data.frame(risk = riskscore[, 1], acture.y)
  ROC <- with(
    timerocDat,
    timeROC::timeROC(
      T = time, delta = status,
      marker = risk, cause = 1,
      weighting = "marginal",
      time = quantile(time, probs = time_prob),
      ROC = TRUE,
      iid = TRUE
    )
  )
  return(ROC)
}

#' Plot Time-Dependent ROC Curves
#'
#' This function generates time-dependent Receiver Operating Characteristic (ROC) curves
#' for evaluating the prognostic accuracy of models based on survival data. It handles
#' both training and testing datasets to plot ROC curves at specific time quantiles.
#' Customizable options for colors and legends are provided via function parameters.
#'
#' @param train.x Matrix or data frame containing the predictor variables used to fit
#'        the model for the training dataset. These variables should correspond to the
#'        same predictors used in model fitting.
#' @param train.y Training dataset outcomes, which must include survival time and
#'        event status (e.g., censoring status). This data should be formatted as a
#'        two-column data frame or matrix where the first column is the survival time
#'        and the second column is the event status.
#' @param test.x Matrix or data frame containing the predictor variables used for
#'        evaluating the model on the testing dataset. Like train.x, this should include
#'        the same type of predictors used in model fitting.
#' @param test.y Testing dataset outcomes, formatted in the same way as train.y, with
#'        survival time and event status.
#' @param model The model object used for predictions. Typically, this would be a model
#'        object created by a survival analysis method compatible with time-dependent
#'        ROC calculations.
#' @param modelname A string representing the name of the model, used for creating
#'        titles or labels in the plots.
#' @param cols Optionally, a vector of colors for plotting. If not provided, colors
#'        are automatically chosen based on the 'palette' parameter.
#' @param palette A character string specifying the color palette to use if 'cols' is
#'        not provided. The default is "jama", which refers to a pre-defined palette.
#'
#' @return A ggplot object representing the ROC curve plot for the provided model at
#'         the specified time quantile. The plot includes both training and testing
#'         datasets across different regularization strengths or model specifications.
#'
#' @examples
#' # Assuming model and data are predefined:
#' PlotTimeROC(train.x, train.y, test.x, test.y, fitted_model, "Cox Model")
#' @export
PlotTimeROC <- function(train.x, train.y, test.x, test.y, model, modelname, cols = NULL, palette = "jama") {
  if (is.null(cols)) {
    cols <- IOBR::palettes(category = "box", palette = palette, show_message = FALSE, show_col = FALSE)
  } else {
    cols <- cols
  }

  newx <- list(train.x, train.x, test.x, test.x)
  s <- list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y <- list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)
  auc <- args %>%
    purrr::pmap(PrognosticAUC, model = model) %>%
    do.call(rbind, .)
  rownames(auc) <- c("Train_lambda.min", "Train_lambda.1se", "Test_lambda.min", "Test_lambda.1se")
  roclist <- args %>% purrr::pmap(CalculateTimeROC, model = model)
  aucs <- round(auc$probs.9, 2)
  legend.name <- paste(c(
    "train_lambda.min", "train_lambda.1se",
    "test_lambda.min", "test_lambda.1se"
  ), "AUC", aucs, sep = " ")
  names(roclist) <- c(
    "train_lambda.min", "train_lambda.1se",
    "test_lambda.min", "test_lambda.1se"
  )

  plotdat <- lapply(roclist, function(z) {
    data.frame(x = z$FP[, 2], y = z$TP[, 2])
  }) %>% plyr::ldply(., .fun = "rbind", .id = "s")
  plotdat$s <- factor(plotdat$s, levels = names(roclist))

  p <- ggplot2::ggplot(plotdat, aes(x = x, y = y)) +
    geom_path(aes(color = s)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("False positive rate") +
    ylab("True positive rate") +
    theme_bw() +
    scale_color_manual(
      values = cols,
      labels = legend.name
    ) +
    ggtitle(paste0(str_replace(modelname, "_", " "), "\nROC at time quantile 0.9")) +
    theme(legend.title = element_blank()) +
    theme(
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.text.x = element_text(face = "plain", angle = 0, hjust = 1, color = "black"),
      axis.text.y = element_text(face = "plain", angle = 0, hjust = 1, color = "black")
    )

  # ggplot2::ggsave(paste0(foldername, "/", modelname, "_ROC.pdf"), plot = p, width = 6, height = 4)
  return(p)
}
