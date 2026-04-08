#' Build Prognostic Models Using LASSO and Ridge Regression
#'
#' @description
#' Prepares data, splits it into training and testing sets, and fits LASSO and Ridge
#' regression models for survival analysis. Evaluates model performance using
#' cross-validation and optionally generates time-dependent ROC curves for visual
#' assessment of predictive accuracy.
#'
#' @param x A matrix or data frame of predictor variables (features).
#' @param y A data frame of survival outcomes with two columns: survival time and event status.
#' @param scale Logical indicating whether to scale predictor variables. Default is `FALSE`.
#' @param seed Integer seed for random number generation to ensure reproducibility.
#'   Default is `123456`.
#' @param train_ratio Numeric proportion of data for training (e.g., 0.7). Default is `0.7`.
#' @param nfold Integer number of folds for cross-validation. Default is `10`.
#' @param plot Logical indicating whether to plot ROC curves. Default is `TRUE`.
#' @param cols Optional vector of colors for ROC curves. If `NULL`, uses default palette.
#' @param palette String specifying color palette. Default is `"jama"`.
#'
#' @return A list containing:
#' \describe{
#'   \item{lasso_result}{Results from LASSO model including coefficients and AUC}
#'   \item{ridge_result}{Results from Ridge model including coefficients and AUC}
#'   \item{train.x}{Training data with sample IDs}
#'   \item{plots}{Named list with \code{lasso} and \code{ridge} time-ROC plots}
#' }
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("glmnet", quietly = TRUE) &&
#'   requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   imvigor210_sig <- load_data("imvigor210_sig")
#'   imvigor210_pdata <- load_data("imvigor210_pdata")
#'   pdata_prog <- data.frame(
#'     ID = imvigor210_pdata$ID,
#'     OS_days = as.numeric(imvigor210_pdata$OS_days),
#'     OS_status = as.numeric(imvigor210_pdata$OS_status)
#'   )
#'   prognostic_result <- PrognosticModel(
#'     x = imvigor210_sig, y = pdata_prog,
#'     scale = TRUE, seed = 123456,
#'     train_ratio = 0.7, nfold = 10, plot = FALSE
#'   )
#' }
#' }
PrognosticModel <- function(x, y, scale = FALSE, seed = 123456, train_ratio = 0.7,
                            nfold = 10, plot = TRUE, palette = "jama", cols = NULL) {
  rlang::check_installed("glmnet")

  x <- as.data.frame(x)
  y <- as.data.frame(y)

  cli::cli_alert_info("Processing data")
  processdat <- ProcessingData(x = x, y = y, scale = scale, type = "survival")
  x_scale <- processdat$x_scale
  y <- processdat$y
  x_ID <- processdat$x_ID

  cli::cli_alert_info("Splitting data into training and validation sets")
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

  cli::cli_alert_info("Running LASSO")
  set.seed(seed)
  lasso_model <- glmnet::cv.glmnet(
    x = train.x, y = as.matrix(train.y),
    family = "cox", alpha = 1, nfolds = nfold
  )
  lasso_result <- PrognosticResult(model = lasso_model, train.x, train.y, test.x, test.y)

  p1 <- PlotTimeROC(
    train.x = train.x, train.y = train.y,
    test.x = test.x, test.y = test.y, model = lasso_model, modelname = "LASSO"
  )
  if (plot) print(p1)

  cli::cli_alert_info("Running RIDGE REGRESSION")
  set.seed(seed)
  ridge_model <- glmnet::cv.glmnet(
    x = train.x, y = as.matrix(train.y),
    family = "cox", alpha = 0, nfolds = nfold
  )
  ridge_result <- PrognosticResult(model = ridge_model, train.x, train.y, test.x, test.y)

  p2 <- PlotTimeROC(
    train.x = train.x, train.y = train.y, cols = cols, palette = palette,
    test.x = test.x, test.y = test.y, model = ridge_model, modelname = "RIDGE"
  )
  if (plot) print(p2)

  cli::cli_alert_success("Model fitting complete")

  list(
    lasso_result = lasso_result, ridge_result = ridge_result, train.x = return.x,
    plots = list(lasso = p1, ridge = p2)
  )
}

#' Compute Prognostic Results for Survival Models
#'
#' @description
#' Computes and compiles prognostic results from a survival model fitted with `glmnet`.
#' Extracts model coefficients at optimal lambda values (`lambda.min` and `lambda.1se`)
#' and calculates time-dependent AUC metrics for both training and testing datasets.
#'
#' @param model A fitted survival model object (e.g., from `glmnet::cv.glmnet`).
#' @param train.x Matrix or data frame of training predictors.
#' @param train.y Training dataset survival outcomes (time and status).
#' @param test.x Matrix or data frame of testing predictors.
#' @param test.y Testing dataset survival outcomes (time and status).
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{The fitted model object}
#'   \item{coefs}{Data frame of coefficients at `lambda.min` and `lambda.1se`}
#'   \item{AUC}{Data frame with AUC values for train/test at both lambda values}
#' }
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE) &&
#'   requireNamespace("survival", quietly = TRUE) &&
#'   requireNamespace("timeROC", quietly = TRUE)) {
#'   library(survival)
#'   set.seed(123)
#'   train_x <- matrix(rnorm(100 * 10), ncol = 10)
#'   train_y <- data.frame(time = rexp(100), status = rbinom(100, 1, 0.5))
#'   test_x <- matrix(rnorm(50 * 10), ncol = 10)
#'   test_y <- data.frame(time = rexp(50), status = rbinom(50, 1, 0.5))
#'   fit <- glmnet::cv.glmnet(train_x, Surv(train_y$time, train_y$status), family = "cox")
#'   results <- PrognosticResult(
#'     model = fit, train.x = train_x, train.y = train_y,
#'     test.x = test_x, test.y = test_y
#'   )
#' }
PrognosticResult <- function(model, train.x, train.y, test.x, test.y) {
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

  datasets <- list(
    list(data = train.x, outcome = train.y, lambda = "lambda.min"),
    list(data = train.x, outcome = train.y, lambda = "lambda.1se"),
    list(data = test.x, outcome = test.y, lambda = "lambda.min"),
    list(data = test.x, outcome = test.y, lambda = "lambda.1se")
  )

  auc_list <- lapply(datasets, function(d) {
    PrognosticAUC(model = model, newx = d$data, s = d$lambda, acture.y = d$outcome)
  })

  auc <- dplyr::bind_rows(auc_list)
  rownames(auc) <- c(
    "Train_lambda.min", "Train_lambda.1se",
    "Test_lambda.min", "Test_lambda.1se"
  )

  list(model = model, coefs = coefs, AUC = auc)
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
#' @param s Lambda value for prediction. Can be numeric or `"lambda.min"`/`"lambda.1se"`.
#' @param acture.y Data frame with `time` and `status` columns.
#'
#' @return A data frame with AUC values at 30th (`probs.3`) and 90th (`probs.9`) percentiles.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE) &&
#'   requireNamespace("survival", quietly = TRUE) &&
#'   requireNamespace("timeROC", quietly = TRUE)) {
#'   library(survival)
#'   set.seed(123)
#'   x <- matrix(rnorm(100 * 5), ncol = 5)
#'   y <- Surv(rexp(100), rbinom(100, 1, 0.5))
#'   fit <- glmnet::cv.glmnet(x, y, family = "cox")
#'   acture_y <- data.frame(time = y[, 1], status = y[, 2])
#'   auc_results <- PrognosticAUC(fit, newx = x, s = "lambda.min", acture.y = acture_y)
#' }
PrognosticAUC <- function(model, newx, s, acture.y) {
  rlang::check_installed("timeROC")
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

  data.frame(probs.3 = ROC$AUC[1], probs.9 = ROC$AUC[2])
}

#' Calculate Time-Dependent ROC Curve
#'
#' @description
#' Computes time-dependent ROC curve for survival models using the `timeROC` package.
#' Evaluates predictive accuracy at a specified time quantile.
#'
#' @param model A fitted survival model object.
#' @param newx A matrix or data frame of new data for prediction.
#' @param s Lambda value for prediction.
#' @param acture.y Data frame with `time` and `status` columns.
#' @param modelname Character string for model identification.
#' @param time_prob Numeric quantile for ROC calculation. Default is `0.9`.
#'
#' @return An object of class `timeROC` containing ROC curve information.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE) &&
#'   requireNamespace("survival", quietly = TRUE) &&
#'   requireNamespace("timeROC", quietly = TRUE)) {
#'   library(survival)
#'   dat <- na.omit(lung[, c("time", "status", "age", "sex", "ph.ecog")])
#'   dat$status <- dat$status - 1
#'   x <- as.matrix(dat[, c("age", "sex", "ph.ecog")])
#'   y <- Surv(dat$time, dat$status)
#'   fit <- glmnet::glmnet(x, y, family = "cox")
#'   actual_outcome <- data.frame(time = dat$time, status = dat$status)
#'   roc_info <- CalculateTimeROC(
#'     model = fit, newx = x, s = 0.01, acture.y = actual_outcome,
#'     modelname = "glmnet Cox Model", time_prob = 0.5
#'   )
#'   print(roc_info$AUC)
#' }
CalculateTimeROC <- function(model, newx, s, acture.y, modelname, time_prob = 0.9) {
  rlang::check_installed("timeROC")
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

  ROC
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
#' @param palette Character string specifying color palette. Default is `"jama"`.
#'
#' @return A `ggplot` object representing the ROC curve plot.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE) &&
#'   requireNamespace("survival", quietly = TRUE) &&
#'   requireNamespace("timeROC", quietly = TRUE)) {
#'   library(survival)
#'   set.seed(123)
#'   train_x <- matrix(rnorm(100 * 5), ncol = 5)
#'   train_y <- data.frame(time = rexp(100), status = rbinom(100, 1, 0.5))
#'   test_x <- matrix(rnorm(50 * 5), ncol = 5)
#'   test_y <- data.frame(time = rexp(50), status = rbinom(50, 1, 0.5))
#'   fit <- glmnet::cv.glmnet(train_x, Surv(train_y$time, train_y$status), family = "cox")
#'   p <- PlotTimeROC(train_x, train_y, test_x, test_y, fit, "Cox Model")
#'   print(p)
#' }
PlotTimeROC <- function(train.x, train.y, test.x, test.y, model, modelname,
                        cols = NULL, palette = "jama") {
  if (is.null(cols)) {
    cols <- palettes(category = "box", palette = palette, show_message = FALSE, show_col = FALSE)
  }

  datasets <- list(
    list(data = train.x, outcome = train.y, lambda = "lambda.min"),
    list(data = train.x, outcome = train.y, lambda = "lambda.1se"),
    list(data = test.x, outcome = test.y, lambda = "lambda.min"),
    list(data = test.x, outcome = test.y, lambda = "lambda.1se")
  )

  auc_list <- lapply(datasets, function(d) {
    PrognosticAUC(model = model, newx = d$data, s = d$lambda, acture.y = d$outcome)
  })

  auc <- dplyr::bind_rows(auc_list)
  rownames(auc) <- c(
    "Train_lambda.min", "Train_lambda.1se",
    "Test_lambda.min", "Test_lambda.1se"
  )

  roclist <- lapply(datasets, function(d) {
    CalculateTimeROC(
      model = model, newx = d$data, s = d$lambda,
      acture.y = d$outcome, modelname = modelname
    )
  })

  aucs <- round(auc$probs.9, 2)
  legend.name <- paste(
    c("train_lambda.min", "train_lambda.1se", "test_lambda.min", "test_lambda.1se"),
    "AUC", aucs,
    sep = " "
  )
  names(roclist) <- c(
    "train_lambda.min", "train_lambda.1se",
    "test_lambda.min", "test_lambda.1se"
  )

  plotdat <- do.call(rbind, lapply(names(roclist), function(nm) {
    data.frame(
      s = nm,
      x = roclist[[nm]]$FP[, 2],
      y = roclist[[nm]]$TP[, 2]
    )
  }))

  plotdat$s <- factor(plotdat$s, levels = names(roclist))

  ggplot2::ggplot(plotdat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$s)) +
    ggplot2::geom_path(linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = cols, labels = legend.name) +
    ggplot2::ggtitle(paste0(stringr::str_replace(modelname, "_", " "), "\nROC at time quantile 0.9")) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0.5),
      axis.text.x = ggplot2::element_text(face = "plain", angle = 0, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(face = "plain", angle = 0, hjust = 1, color = "black")
    )
}
