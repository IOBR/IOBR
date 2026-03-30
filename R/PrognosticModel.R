#' Build Prognostic Models Using LASSO and Ridge Regression
#'
#' @description
#' Prepares data, splits it into training and testing sets, and fits LASSO and Ridge regression models
#' for survival analysis. Evaluates model performance using cross-validation and optionally generates
#' time-dependent ROC curves for visual assessment of predictive accuracy.
#'
#' @param x A matrix or data frame of predictor variables (features).
#' @param y A data frame of survival outcomes with two columns: survival time and event status.
#' @param scale Logical indicating whether to scale predictor variables. Default is \code{FALSE}.
#' @param seed Integer seed for random number generation to ensure reproducibility. Default is \code{123456}.
#' @param train_ratio Numeric proportion of data for training (e.g., 0.7). Default is \code{0.7}.
#' @param nfold Integer number of folds for cross-validation. Default is \code{10}.
#' @param plot Logical indicating whether to plot ROC curves. Default is \code{TRUE}.
#' @param cols Optional vector of colors for ROC curves. If \code{NULL}, uses default palette.
#' @param palette String specifying color palette. Default is \code{"jama"}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{lasso_result}: Results from LASSO model including coefficients and AUC
#'   \item \code{ridge_result}: Results from Ridge model including coefficients and AUC
#'   \item \code{train.x}: Training data with sample IDs
#' }
#' @author Dongqiang Zeng
#' @export
#' @importFrom glmnet cv.glmnet
#' @import dplyr
#' @import ggplot2
#' @examples
#' \dontrun{
#' imvigor210_sig <- load_data("imvigor210_sig")
#' imvigor210_pdata <- load_data("imvigor210_pdata")
#' pdata_prog <- imvigor210_pdata %>%
#'   dplyr::select(ID, OS_days, OS_status) %>%
#'   dplyr::mutate(OS_days = as.numeric(OS_days), OS_status = as.numeric(OS_status))
#' prognostic_result <- PrognosticModel(
#'   x = imvigor210_sig, y = pdata_prog,
#'   scale = TRUE, seed = 123456,
#'   train_ratio = 0.7, nfold = 10, plot = TRUE
#' )
#' }
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
  if (plot) {
    p <- p1 + p2
    print(p)
  } else {
    p <- NULL
  }
  return(list(lasso_result = lasso_result, ridge_result = ridge_result, train.x = return.x))
}

#' Compute Prognostic Results for Survival Models
#'
#' @description
#' Computes and compiles prognostic results from a survival model fitted with \code{glmnet}.
#' Extracts model coefficients at optimal lambda values (\code{lambda.min} and \code{lambda.1se})
#' and calculates time-dependent AUC metrics for both training and testing datasets.
#'
#' @param model A fitted survival model object (e.g., from \code{glmnet::cv.glmnet}).
#' @param train.x Matrix or data frame of training predictors.
#' @param train.y Training dataset survival outcomes (time and status).
#' @param test.x Matrix or data frame of testing predictors.
#' @param test.y Testing dataset survival outcomes (time and status).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{model}: The fitted model object
#'   \item \code{coefs}: Data frame of coefficients at \code{lambda.min} and \code{lambda.1se}
#'   \item \code{AUC}: Data frame with AUC values for train/test at both lambda values
#' }
#'
#' @author Dongqiang Zeng
#' @export
#' @importFrom stats coef
#' @importFrom purrr pmap
#' @examples
#' \dontrun{
#' # Example with a fitted glmnet Cox model
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(123)
#'   train_data <- list(
#'     x = matrix(rnorm(100 * 10), ncol = 10),
#'     y = survival::Surv(rexp(100), rbinom(100, 1, 0.5))
#'   )
#'   test_data <- list(
#'     x = matrix(rnorm(50 * 10), ncol = 10),
#'     y = survival::Surv(rexp(50), rbinom(50, 1, 0.5))
#'   )
#'   fit <- glmnet::cv.glmnet(train_data$x, train_data$y, family = "cox")
#'   results <- PrognosticResult(
#'     model = fit, train.x = train_data$x, train.y = train_data$y,
#'     test.x = test_data$x, test.y = test_data$y
#'   )
#' }
#' }
#' @export
PrognosticResult <- function(model, train.x, train.y, test.x, test.y) {
  # Extract coefficients
  coefs <- cbind(
    stats::coef(model, s = "lambda.min"),
    stats::coef(model, s = "lambda.1se")
  )
  coefs <- data.frame(
    feature = rownames(coefs),
    lambda.min = coefs[, 1],
    lambda.1se = coefs[, 2],
    row.names = NULL
  )

  # Prepare AUC calculation parameters
  args <- list(
    newx = list(train.x, train.x, test.x, test.x),
    s = list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se"),
    acture.y = list(train.y, train.y, test.y, test.y)
  )

  # Explicitly pass parameters to avoid pmap misalignment
  auc_list <- purrr::pmap(
    args,
    function(newx, s, acture.y) {
      PrognosticAUC(model = model, newx = newx, s = s, acture.y = acture.y)
    }
  )

  auc <- do.call(rbind, auc_list)
  rownames(auc) <- c(
    "Train_lambda.min", "Train_lambda.1se",
    "Test_lambda.min", "Test_lambda.1se"
  )

  resultreturn <- list(model = model, coefs = coefs, AUC = auc)
  return(resultreturn)
}


#' Calculate Time-Dependent AUC for Survival Models
#'
#' @description
#' Evaluates prognostic ability of a survival model by calculating time-dependent AUC
#' at the 30th and 90th percentiles of survival time. These thresholds assess
#' short-term and long-term predictive accuracy.
#'
#' @param model A fitted survival model object capable of generating risk scores.
#' @param newx A matrix or data frame of new data for prediction.
#' @param s Lambda value for prediction. Can be numeric or \code{"lambda.min"}/\code{"lambda.1se"}.
#' @param acture.y Data frame with \code{time} and \code{status} columns.
#'
#' @return A data frame with AUC values at 30th (\code{probs.3}) and 90th (\code{probs.9}) percentiles.
#' @author Dongqiang Zeng
#' @export
#' @importFrom stats predict
#' @importFrom timeROC timeROC
#' @examples
#' \dontrun{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(123)
#'   x <- matrix(rnorm(100 * 5), ncol = 5)
#'   y <- survival::Surv(rexp(100), rbinom(100, 1, 0.5))
#'   fit <- glmnet::cv.glmnet(x, y, family = "cox")
#'   acture_y <- data.frame(time = y[, 1], status = y[, 2])
#'   auc_results <- PrognosticAUC(fit, newx = x, s = "lambda.min", acture.y = acture_y)
#' }
#' }
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


#' Calculate Time-Dependent ROC Curve
#'
#' @description
#' Computes time-dependent ROC curve for survival models using the \code{timeROC} package.
#' Evaluates predictive accuracy at a specified time quantile.
#'
#' @param model A fitted survival model object.
#' @param newx A matrix or data frame of new data for prediction.
#' @param s Lambda value for prediction.
#' @param acture.y Data frame with \code{time} and \code{status} columns.
#' @param modelname Character string for model identification.
#' @param time_prob Numeric quantile for ROC calculation (default: \code{0.9}).
#'
#' @return An object of class \code{timeROC} containing ROC curve information.
#' @author Dongqiang Zeng
#' @export
#' @importFrom stats predict
#' @importFrom timeROC timeROC
#' @examples
#' \dontrun{
#' if (requireNamespace("glmnet", quietly = TRUE) &&
#'     requireNamespace("survival", quietly = TRUE)) {
#'   # Use lung dataset and remove missing values
#'   dat <- na.omit(survival::lung[, c("time", "status", "age", "sex", "ph.ecog")])
#'   # Convert status to 0/1 for glmnet Cox model
#'   dat$status <- dat$status - 1
#'   x <- as.matrix(dat[, c("age", "sex", "ph.ecog")])
#'   y <- survival::Surv(dat$time, dat$status)
#'   fit <- glmnet::glmnet(x, y, family = "cox")
#'   actual_outcome <- data.frame(time = dat$time, status = dat$status)
#'   roc_info <- CalculateTimeROC(
#'     model = fit, newx = x, s = 0.01, acture.y = actual_outcome,
#'     modelname = "glmnet Cox Model", time_prob = 0.5
#'   )
#'   print(roc_info$AUC)
#' }
#' }
#' @export
CalculateTimeROC <- function(model, newx, s, acture.y, modelname, time_prob = 0.9) {
  riskscore <- stats::predict(model, newx = newx, s = s)
  timerocDat <- data.frame(risk = riskscore[, 1], acture.y)
  # Ensure survival package is available
  if (!exists("Surv", mode = "function")) {
    stop("Surv function not available. Please ensure survival package is loaded.")
  }
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
#' @description
#' Generates time-dependent ROC curves for evaluating prognostic accuracy of survival models.
#' Plots training and testing ROC curves at the 90th percentile survival time.
#'
#' @param train.x Matrix or data frame of training predictors.
#' @param train.y Training survival outcomes (time and status).
#' @param test.x Matrix or data frame of testing predictors.
#' @param test.y Testing survival outcomes (time and status).
#' @param model Fitted survival model object.
#' @param modelname Character string for model identification.
#' @param cols Optional vector of colors for plotting.
#' @param palette Character string specifying color palette (default: \code{"jama"}).
#'
#' @return A \code{ggplot} object representing the ROC curve plot.
#' @author Dongqiang Zeng
#' @export
#' @examples
#' \dontrun{
#' # Example usage with fitted model
#' PlotTimeROC(train.x, train.y, test.x, test.y, fitted_model, "Cox Model")
#' }
#' @export
PlotTimeROC <- function(train.x, train.y, test.x, test.y, model, modelname, cols = NULL, palette = "jama") {
  if (is.null(cols)) {
    cols <- palettes(category = "box", palette = palette, show_message = FALSE, show_col = FALSE)
  } else {
    cols <- cols
  }

  args <- list(
    newx = list(train.x, train.x, test.x, test.x),
    s = list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se"),
    acture.y = list(train.y, train.y, test.y, test.y)
  )

  auc_list <- purrr::pmap(
    args,
    function(newx, s, acture.y) {
      PrognosticAUC(model = model, newx = newx, s = s, acture.y = acture.y)
    }
  )
  auc <- do.call(rbind, auc_list)
  rownames(auc) <- c(
    "Train_lambda.min", "Train_lambda.1se",
    "Test_lambda.min", "Test_lambda.1se"
  )

  roclist <- purrr::pmap(
    args,
    function(newx, s, acture.y) {
      CalculateTimeROC(
        model = model,
        newx = newx,
        s = s,
        acture.y = acture.y,
        modelname = modelname
      )
    }
  )

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
    ggtitle(paste0(stringr::str_replace(modelname, "_", " "), "\nROC at time quantile 0.9")) +
    theme(legend.title = element_blank()) +
    theme(
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.text.x = element_text(face = "plain", angle = 0, hjust = 1, color = "black"),
      axis.text.y = element_text(face = "plain", angle = 0, hjust = 1, color = "black")
    )

  # ggplot2::ggsave(paste0(foldername, "/", modelname, "_ROC.pdf"), plot = p, width = 6, height = 4)
  return(p)
}
