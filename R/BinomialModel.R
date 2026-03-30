#' Binomial Model Construction
#'
#' @description
#' Constructs and evaluates binomial logistic regression models using Lasso and Ridge
#' regularization. Processes input data, scales features if specified, splits data into
#' training/testing sets, and fits both Lasso and Ridge models. Optionally generates
#' AUC plots for model evaluation.
#'
#' @param x A data frame containing sample ID and features. First column must be sample ID.
#' @param y A data frame where first column is sample ID and second column is outcome
#'   (numeric or factor).
#' @param seed Integer for random seed. Default is `123456`.
#' @param scale Logical indicating whether to scale features. Default is `TRUE`.
#' @param train_ratio Numeric between 0 and 1 for training proportion. Default is `0.7`.
#' @param nfold Integer for cross-validation folds. Default is `10`.
#' @param plot Logical indicating whether to generate AUC plots. Default is `TRUE`.
#' @param cols Optional color vector for ROC curves. Default is `NULL`.
#' @param palette Character string for color palette. Default is `"jama"`.
#'
#' @return List containing:
#' \describe{
#'   \item{lasso_result}{Lasso model results}
#'   \item{ridge_result}{Ridge model results}
#'   \item{train.x}{Training data with IDs}
#' }
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' imvigor210_sig <- load_data("imvigor210_sig")
#' imvigor210_pdata <- load_data("imvigor210_pdata")
#' pdata_group <- imvigor210_pdata[imvigor210_pdata$BOR_binary != "NA", c("ID", "BOR_binary")]
#' pdata_group$BOR_binary <- factor(ifelse(pdata_group$BOR_binary == "R", 1, 0))
#' result <- BinomialModel(
#'   x = imvigor210_sig, y = pdata_group,
#'   seed = 123456, scale = TRUE, train_ratio = 0.7, nfold = 10, plot = FALSE
#' )
#' }
BinomialModel <- function(x, y, seed = 123456, scale = TRUE, train_ratio = 0.7,
                         nfold = 10, plot = TRUE, palette = "jama", cols = NULL) {

  rlang::check_installed("glmnet")

  x <- as.data.frame(x)
  y <- as.data.frame(y)

  cli::cli_alert_info("Processing data")
  processdat <- ProcessingData(x = x, y = y, scale = scale, type = "binomial")
  x_scale <- processdat$x_scale
  y <- processdat$y
  x_ID <- processdat$x_ID

  cli::cli_alert_info("Splitting data into training and test sets")
  train_test <- SplitTrainTest(
    x = x_scale, y = y, train_ratio = train_ratio,
    type = "binomial", seed = seed
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
    x = train.x, y = train.y, family = "binomial",
    type.measure = "class", alpha = 1, nfolds = nfold
  )
  lasso_result <- RegressionResult(
    train.x = train.x, train.y = train.y,
    test.x = test.x, test.y = test.y, model = lasso_model
  )

  if (plot) {
    p1 <- PlotAUC(
      train.x = train.x, train.y = train.y,
      test.x = test.x, test.y = test.y, model = lasso_model,
      cols = cols, palette = palette, modelname = "lasso_model"
    )
    print(p1)
  }

  cli::cli_alert_info("Running RIDGE REGRESSION")
  set.seed(seed)
  ridge_model <- glmnet::cv.glmnet(
    x = train.x, y = train.y, family = "binomial",
    type.measure = "class", alpha = 0, nfolds = nfold
  )
  ridge_result <- RegressionResult(
    train.x = train.x, train.y = train.y,
    test.x = test.x, test.y = test.y, model = ridge_model
  )

  if (plot) {
    p2 <- PlotAUC(
      train.x = train.x, train.y = train.y,
      test.x = test.x, test.y = test.y, model = ridge_model,
      cols = cols, palette = palette, modelname = "ridge_model"
    )
    print(p2)
  }

  cli::cli_alert_success("Model fitting complete")

  list(lasso_result = lasso_result, ridge_result = ridge_result, train.x = return.x)
}

#' Process Data for Model Construction
#'
#' @description
#' Preprocesses data for binomial or survival analysis. Aligns and filters data
#' based on sample IDs, optionally scales data, and ensures appropriate data types.
#'
#' @param x Data frame of predictors with first column as IDs.
#' @param y Data frame of outcomes with first column as IDs. For survival,
#'   expects two additional columns for time and status.
#' @param scale Logical indicating whether to scale predictors.
#' @param type Character string: `"binomial"` or `"survival"`.
#'
#' @return List containing:
#' \describe{
#'   \item{x_scale}{Processed predictor matrix}
#'   \item{y}{Processed outcome variable}
#'   \item{x_ID}{Sample IDs}
#' }
#'
#' @export
#' @examples
#' \donttest{
#' x <- data.frame(ID = 1:10, predictor1 = rnorm(10), predictor2 = rnorm(10))
#' y <- data.frame(ID = 1:10, outcome = sample(c(0, 1), 10, replace = TRUE))
#' result <- ProcessingData(x, y, scale = TRUE, type = "binomial")
#' }
ProcessingData <- function(x, y, scale, type = c("binomial", "survival")) {

  type <- rlang::arg_match(type)

  colnames(x)[1] <- "ID"
  colnames(y)[1] <- "ID"
  x$ID <- as.character(x$ID)
  y$ID <- as.character(y$ID)

  if (type == "survival") {
    colnames(y) <- c("ID", "time", "status")
    y <- dplyr::filter(y, .data$time > 0)
  }

  samples <- intersect(x$ID, y$ID)
  if (length(samples) == 0) {
    cli::cli_abort("No matching sample IDs found between x and y")
  }

  x <- x[match(samples, x$ID), ]
  y <- y[match(samples, y$ID), ]

  if (type == "binomial") {
    colnames(y) <- c("ID", "Group")
    y <- dplyr::pull(y, .data$Group)
    if (!is.factor(y)) {
      cli::cli_warn("Converting outcome to factor")
      y <- as.factor(y)
    }
  } else if (type == "survival") {
    y <- y[, c("time", "status")]
  }

  x_scale <- if (scale) {
    scale(x[, -1], center = TRUE, scale = TRUE)
  } else {
    x[, -1]
  }

  x_ID <- x[, "ID"]
  na_cols <- which(as.numeric(apply(x_scale, 2, function(z) sum(is.na(z)))) != 0)
  if (length(na_cols) > 0) {
    x_scale <- x_scale[, -na_cols]
  }

  list(x_scale = x_scale, y = y, x_ID = x_ID)
}

#' Regression Result Computation
#'
#' @description
#' Computes regression results with coefficients at lambda.min and lambda.1se,
#' and evaluates AUC for binomial outcomes.
#'
#' @param train.x Training predictors matrix.
#' @param train.y Training outcomes.
#' @param test.x Testing predictors matrix.
#' @param test.y Testing outcomes.
#' @param model Fitted model object.
#'
#' @return List containing model, coefficients, and AUC values.
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(123)
#'   train_data <- matrix(rnorm(100 * 10), ncol = 10)
#'   train_outcome <- rbinom(100, 1, 0.5)
#'   test_data <- matrix(rnorm(50 * 10), ncol = 10)
#'   test_outcome <- rbinom(50, 1, 0.5)
#'   fitted_model <- glmnet::cv.glmnet(train_data, train_outcome, family = "binomial", nfolds = 5)
#'   results <- RegressionResult(
#'     train.x = train_data, train.y = train_outcome,
#'     test.x = test_data, test.y = test_outcome, model = fitted_model
#'   )
#' }
#' }
RegressionResult <- function(train.x, train.y, test.x, test.y, model) {
  coefs <- cbind(
    stats::coef(model, s = model$lambda.min),
    stats::coef(model, s = model$lambda.1se)
  )
  coefs <- data.frame(
    feature = rownames(coefs),
    lambda.min = coefs[, 1],
    lambda.1se = coefs[, 2]
  )

  args <- list(
    newx = list(train.x, train.x, test.x, test.x),
    s = list(model$lambda.min, model$lambda.1se, model$lambda.min, model$lambda.1se),
    acture.y = list(train.y, train.y, test.y, test.y)
  )

  aucs <- purrr::pmap_dbl(args, BinomialAUC, model = model)
  AUC <- matrix(aucs, ncol = 2, byrow = TRUE,
                dimnames = list(c("train", "test"), c("lambda.min", "lambda.1se")))

  list(model = model, coefs = coefs, AUC = AUC)
}

#' Elastic Net Model Fitting
#'
#' @description
#' Fits elastic net model with cross-validation to find optimal alpha and lambda.
#'
#' @param train.x Training predictors matrix.
#' @param train.y Training binary outcomes.
#' @param lambdamax Maximum lambda value.
#' @param nfold Number of CV folds. Default is `10`.
#'
#' @return List with `chose_alpha` and `chose_lambda`.
#'
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(123)
#'   train_data <- matrix(rnorm(50 * 5), ncol = 5)
#'   train_outcome <- rbinom(50, 1, 0.5)
#'   result <- Enet(train.x = train_data, train.y = train_outcome, lambdamax = 1, nfold = 5)
#' }
#' }
Enet <- function(train.x, train.y, lambdamax, nfold = 10) {
  train.x <- as.matrix(train.x)
  train.y <- as.numeric(train.y)

  if (nrow(train.x) != length(train.y)) {
    cli::cli_abort("Number of rows in train.x must equal length of train.y")
  }
  if (lambdamax <= 0) {
    cli::cli_abort("lambdamax must be greater than 0")
  }
  if (nfold < 2) {
    cli::cli_abort("nfold must be at least 2")
  }

  alpha.grid <- seq(0, 1, by = 0.2)
  lambda.grid <- seq(1e-4, lambdamax, length.out = 10)

  res <- lapply(alpha.grid, function(a) {
    fit <- glmnet::cv.glmnet(
      x = train.x, y = train.y, family = "binomial",
      alpha = a, lambda = lambda.grid, nfolds = nfold, type.measure = "class"
    )
    data.frame(alpha = a, lambda = fit$lambda.min, cv_error = min(fit$cvm))
  })

  res <- do.call(rbind, res)
  best <- res[which.min(res$cv_error), ]

  list(chose_alpha = best$alpha, chose_lambda = best$lambda)
}

#' Calculate AUC for Binomial Model
#'
#' @description
#' Computes AUC for model predictions.
#'
#' @param model Fitted model object.
#' @param newx New data for prediction.
#' @param s Lambda value for prediction.
#' @param acture.y Actual binary outcomes.
#'
#' @return Numeric AUC value.
#'
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("glmnet", quietly = TRUE) && requireNamespace("ROCR", quietly = TRUE)) {
#'   fitted_model <- glmnet::glmnet(matrix(rnorm(100), ncol = 2), rbinom(50, 1, 0.5), family = "binomial")
#'   auc_value <- BinomialAUC(fitted_model, matrix(rnorm(20), ncol = 2), "lambda.min", rbinom(10, 1, 0.5))
#' }
#' }
BinomialAUC <- function(model, newx, s, acture.y) {
  rlang::check_installed("ROCR")
  prob <- stats::predict(model, newx = newx, s = s, type = "response")
  pred <- ROCR::prediction(prob, acture.y)
  as.numeric(ROCR::performance(pred, "auc")@y.values)
}

#' Plot AUC ROC Curves
#'
#' @description
#' Generates ROC curves for model evaluation.
#'
#' @param train.x Training predictors.
#' @param train.y Training outcomes.
#' @param test.x Testing predictors.
#' @param test.y Testing outcomes.
#' @param model Fitted model.
#' @param modelname Model name for title.
#' @param cols Optional color vector.
#' @param palette Color palette name.
#'
#' @return ggplot object of ROC curve.
#'
#' @export
#' @examples
#' \dontrun{
#' PlotAUC(train.x, train.y, test.x, test.y, fitted_model, "MyModel")
#' }
PlotAUC <- function(train.x, train.y, test.x, test.y, model, modelname,
                   cols = NULL, palette = "jama") {

  cols <- cols %||% palettes(category = "box", palette = palette,
                            show_message = FALSE, show_col = FALSE)

  args <- list(
    newx = list(train.x, train.x, test.x, test.x),
    s = list(model$lambda.min, model$lambda.1se, model$lambda.min, model$lambda.1se),
    acture.y = list(train.y, train.y, test.y, test.y)
  )

  pref <- purrr::pmap(args, CalculatePref, model = model)
  aucs <- round(purrr::pmap_dbl(args, BinomialAUC, model = model), 2)

  legend.name <- paste(
    c("train_lambda.min", "train_lambda.1se", "test_lambda.min", "test_lambda.1se"),
    "AUC", aucs, sep = " "
  )
  names(pref) <- c("train_lambda.min", "train_lambda.1se", "test_lambda.min", "test_lambda.1se")

  plotdat <- purrr::map_dfr(pref, function(z) {
    data.frame(x = z@x.values[[1]], y = z@y.values[[1]])
  }, .id = "s")

  plotdat$s <- factor(plotdat$s, levels = names(pref))

  ggplot2::ggplot(plotdat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$s)) +
    ggplot2::geom_path(linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = cols, labels = legend.name) +
    ggplot2::ggtitle(stringr::str_replace(modelname, "_", " ")) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = ggplot2::rel(2), hjust = 0.5),
      axis.text.x = ggplot2::element_text(face = "plain", angle = 0, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(face = "plain", angle = 0, hjust = 1, color = "black")
    )
}

#' Calculate Performance Metrics
#'
#' @description
#' Computes TPR and FPR for ROC analysis.
#'
#' @param model Fitted model.
#' @param newx New data for prediction.
#' @param s Lambda value.
#' @param acture.y Actual outcomes.
#'
#' @return ROCR performance object.
#'
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("glmnet", quietly = TRUE) && requireNamespace("ROCR", quietly = TRUE)) {
#'   fitted_model <- glmnet::cv.glmnet(matrix(rnorm(100), ncol = 2), rbinom(50, 1, 0.5), nfolds = 3)
#'   perf <- CalculatePref(fitted_model, matrix(rnorm(20), ncol = 2), "lambda.min", rbinom(10, 1, 0.5))
#' }
#' }
CalculatePref <- function(model, newx, s, acture.y) {
  rlang::check_installed("ROCR")
  prob <- stats::predict(model, newx = newx, s = s, type = "response")
  pred <- ROCR::prediction(prob, acture.y)
  ROCR::performance(pred, "tpr", "fpr")
}

#' Split Data into Training and Testing Sets
#'
#' @description
#' Divides dataset into training and testing sets.
#'
#' @param x Predictor matrix or data frame.
#' @param y Outcome vector or matrix.
#' @param train_ratio Proportion for training (0-1).
#' @param type `"binomial"` or `"survival"`.
#' @param seed Random seed.
#'
#' @return List with train.x, train.y, test.x, test.y, train_sample.
#'
#' @export
#' @examples
#' \donttest{
#' data_matrix <- matrix(rnorm(200), ncol = 2)
#' outcome_vector <- rbinom(100, 1, 0.5)
#' split_data <- SplitTrainTest(data_matrix, outcome_vector, train_ratio = 0.7, type = "binomial", seed = 123)
#' }
SplitTrainTest <- function(x, y, train_ratio, type = c("binomial", "survival"), seed) {

  type <- rlang::arg_match(type)

  sizes <- round(nrow(x) * train_ratio)
  set.seed(seed)
  train_sample <- sample(seq_len(nrow(x)), size = sizes, replace = FALSE)
  test_sample <- setdiff(seq_len(nrow(x)), train_sample)

  train.x <- x[train_sample, , drop = FALSE]
  test.x <- x[test_sample, , drop = FALSE]

  if (type == "binomial") {
    train.y <- y[train_sample]
    test.y <- y[test_sample]
  } else if (type == "survival") {
    train.y <- y[train_sample, , drop = FALSE]
    test.y <- y[test_sample, , drop = FALSE]
  }

  list(train.x = train.x, train.y = train.y, test.x = test.x, test.y = test.y,
       train_sample = train_sample)
}
