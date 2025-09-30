#' Binomial Model Construction
#'
#' Constructs and evaluates binomial logistic regression models using Lasso and Ridge regularization.
#' This function processes the input data, scales features if specified, splits the data into training
#' and testing sets based on a provided ratio, and then fits both Lasso and Ridge models. It optionally
#' generates AUC plots for model evaluation. The results, along with the training dataset, are returned.
#'
#' @param x A data.frame containing the sample ID and features, the first column must be the sample ID.
#' @param y A data.frame where the first column is the sample ID and the second column is the outcome
#'          of each sample, which can be either numeric or a factor vector.
#' @param seed An integer used to set the seed for random operations, default is 123456.
#' @param scale A logical indicating whether the feature data `x` should be scaled, default is TRUE.
#' @param train_ratio A numeric value between 0 and 1 specifying the proportion of the dataset to be
#'                    used for training, e.g., 0.7 means 70 percent of the data is used for training.
#' @param nfold The number of folds to use for cross-validation in model fitting, default is 10.
#' @param plot A logical indicating whether to generate and display AUC plots, default is TRUE.
#' @param cols Optional, a vector of colors for the ROC curves. If NULL, default color palettes are applied based on the 'palette' parameter.
#' @param palette Optional, a string specifying the color palette. Default is "jama", which is applied if 'cols' is NULL.
#'
#' @return A list containing the results of the Lasso and Ridge models, and the input training data.
#'         The list includes elements 'lasso_result', 'ridge_result', and 'train.x'.
#'
#' @examples
#' data("imvigor210_sig", package = "IOBR")
#' data("imvigor210_pdata", package = "IOBR")
#' pdata_group <- imvigor210_pdata[imvigor210_pdata$BOR_binary != "NA", c("ID", "BOR_binary")]
#' pdata_group$BOR_binary <- ifelse(pdata_group$BOR_binary == "R", 1, 0)
#' BinomialModel(x = imvigor210_sig, y = pdata_group, seed = 123456, scale = TRUE, train_ratio = 0.7, nfold = 10, plot = T)
#' @export
BinomialModel <- function(x, y, seed = 123456, scale = TRUE, train_ratio = 0.7, nfold = 10, plot = TRUE, palette = "jama", cols = NULL) {
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  print(message(paste0("\n", ">>> Processing data")))
  processdat <- ProcessingData(x = x, y = y, scale = scale, type = "binomial")
  x_scale <- processdat$x_scale
  y <- processdat$y
  x_ID <- processdat$x_ID
  print(message(paste0("\n", ">>> Spliting data into train and test data")))
  train_test <- SplitTrainTest(
    x = x_scale, y = y, train_ratio = train_ratio, type = "binomial",
    seed = seed
  )
  train.x <- train_test$train.x
  train.y <- train_test$train.y
  test.x <- train_test$test.x
  test.y <- train_test$test.y
  train_sample <- train_test$train_sample
  return.x <- data.frame(ID = x_ID[train_sample], train.x)

  print(message(paste0("\n", ">>> Running ", "LASSO")))
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
      test.x = test.x, test.y = test.y, model = lasso_model, cols = cols, palette = palette,
      modelname = "lasso_model"
    )
    print(p1)
  }

  print(message(paste0("\n", ">>> Running ", "RIDGE REGRESSION")))

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
      cols = cols, palette = palette,
      modelname = "ridge_model"
    )
    print(p2)
  }
  print(message(paste0("\n", ">>> Running ", "Elastic Network.")))

  message(paste0("\n", ">>> Done !"))
  return(list(
    lasso_result = lasso_result, ridge_result = ridge_result,
    train.x = return.x
  ))
}


#' Processing Data for Model construction
#'
#' This function preprocesses the data for binomial or survival analysis. It aligns and filters the data based on sample IDs,
#' optionally scales the data, and makes sure the data types are appropriate for further modeling. It also handles survival data by filtering out
#' entries where the time is non-positive and adjusts the data structure based on the type of analysis.
#'
#' @param x A data frame containing predictors with the first column being IDs.
#' @param y A data frame containing the outcome variable with the first column being IDs. For survival analysis, it expects two additional columns for time and status.
#' @param scale Logical, indicating whether the predictor variables should be scaled. When TRUE, predictors are centered and scaled.
#' @param type Character string specifying the type of analysis. Possible values are "binomial" for binomial outcomes where y is a factor indicating group membership, or "survival" where y includes survival time and event status.
#'
#' @return A list containing three elements:
#'   - `x_scale`: The processed predictor matrix with optional scaling applied.
#'   - `y`: The outcome variable, processed according to the specified analysis type.
#'   - `x_ID`: The IDs of the samples included in the analysis.
#'
#' @examples
#' # For binomial analysis
#' x <- data.frame(ID = 1:10, predictor1 = rnorm(10), predictor2 = rnorm(10))
#' y <- data.frame(ID = 1:10, outcome = sample(c(0, 1), 10, replace = TRUE))
#' result <- ProcessingData(x, y, scale = TRUE, type = "binomial")
#'
#' # For survival analysis
#' y <- data.frame(ID = 1:10, time = runif(10, 0, 100), status = sample(c(0, 1), 10, replace = TRUE))
#' result <- ProcessingData(x, y, scale = FALSE, type = "survival")
#' @export
ProcessingData <- function(x, y, scale, type = "binomial") {
  require(dplyr)

  colnames(x)[1] <- "ID"
  colnames(y)[1] <- "ID"
  x$ID <- as.character(x$ID)
  y$ID <- as.character(y$ID)

  if (type == "survival") {
    colnames(y) <- c("ID", "time", "status")
    y <- dplyr::filter(y, time > 0)
  }

  samples <- intersect(x$ID, y$ID)
  if (length(samples) == 0) {
    stop("No matching sample ID found between input matrices x and y.")
  }

  x <- x[match(samples, x$ID), ]
  y <- y[match(samples, y$ID), ]

  if (type == "binomial") {
    colnames(y) <- c("ID", "Group")
    y <- dplyr::pull(y, Group)
    if (!is.factor(y)) {
      warning("Outcome is not a factor, transforming it into a factor vector.")
      y <- as.factor(y)
    }
  }

  if (type == "survival") {
    y <- y[, c("time", "status")]
  }

  if (scale) {
    x_scale <- scale(x[, -1], center = TRUE, scale = TRUE)
  } else {
    x_scale <- x[, -1]
  }

  x_ID <- x[, "ID"]
  ValueNA <- which(as.numeric(apply(x_scale, 2, function(z) sum(is.na(z)))) != 0)
  if (length(ValueNA) > 0) {
    x_scale <- x_scale[, -ValueNA]
  }

  return(list(x_scale = x_scale, y = y, x_ID = x_ID))
}


#' Regression Result Computation
#'
#' Computes the regression results for given training and testing datasets using a specified model.
#' It calculates the coefficients at specified regularization strengths (lambda.min and lambda.1se) and
#' evaluates the model's performance using Area Under the Curve (AUC) for binomial outcomes.
#'
#' @param train.x A matrix or data frame of training predictors used to fit the model.
#' @param train.y Training dataset outcomes, expected to be a binary vector indicating the presence or absence of an event.
#' @param test.x A matrix or data frame of testing predictors used for evaluating the model.
#' @param test.y Testing dataset outcomes, expected to be similar in format to train.y.
#' @param model A model object from which coefficients will be extracted and used for evaluation.
#'
#' @return A list containing the model object, a data frame of coefficients for different lambda values,
#'         and a matrix of AUC values for both the training and testing sets, calculated at both lambda.min and lambda.1se.
#'
#' @examples
#' # Assuming that `model` is already fitted using glmnet or a similar package:
#' train_data <- matrix(rnorm(100 * 10), ncol = 10)
#' train_outcome <- rbinom(100, 1, 0.5)
#' test_data <- matrix(rnorm(100 * 10), ncol = 10)
#' test_outcome <- rbinom(100, 1, 0.5)
#' fitted_model <- glmnet(train_data, train_outcome, family = "binomial")
#' results <- RegressionResult(
#'   train.x = train_data, train.y = train_outcome,
#'   test.x = test_data, test.y = test_outcome,
#'   model = fitted_model
#' )
#' @export
RegressionResult <- function(train.x, train.y, test.x, test.y, model) {
  coefs <- cbind(coef(model, s = "lambda.min"), coef(model, s = "lambda.1se"))
  coefs <- data.frame(feature = rownames(coefs), lambda.min = coefs[, 1], lambda.1se = coefs[, 2])
  newx <- list(train.x, train.x, test.x, test.x)
  s <- list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y <- list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)
  AUC <- args %>%
    purrr::pmap_dbl(BinomialAUC, model = model) %>%
    matrix(.,
      ncol = 2, byrow = T,
      dimnames = list(c("train", "test"), c("lambda.min", "lambda.1se"))
    )
  resultreturn <- list(
    model = model, coefs = coefs,
    AUC = AUC
  )
}

#' Elastic Net Model Fitting
#'
#' Fits an elastic net model using the `caret` package, allowing for both L1 (Lasso) and L2 (Ridge) regularization.
#' The function uses cross-validation to identify the optimal alpha and lambda values that maximize accuracy.
#' Alpha ranges from 0 (Ridge) to 1 (Lasso), with lambda controlling the strength of the regularization.
#'
#' @param train.x A matrix or data frame containing training predictors.
#' @param train.y A numeric vector or factor representing the binary outcome for each sample in the training set.
#' @param lambdamax The maximum value of lambda to consider in the regularization path.
#' @param nfold The number of folds to use for cross-validation, which also determines the number of repeats for repeated CV.
#'
#' @return A list containing the optimal values of alpha and lambda chosen based on cross-validation.
#'
#' @examples
#' # Assuming 'train_data' and 'train_outcome' are already defined:
#' train_data <- matrix(rnorm(100 * 10), ncol = 10)
#' train_outcome <- rbinom(100, 1, 0.5)
#' optimal_parameters <- Enet(
#'   train.x = train_data, train.y = train_outcome,
#'   lambdamax = 1, nfold = 10
#' )
#' print(optimal_parameters)
#' @export
Enet <- function(train.x, train.y, lambdamax, nfold = nfold) {
  grid <- expand.grid(.alpha = seq(0, 1, by = .2), .lambda = seq(0, lambdamax, length.out = 10))
  fitControl <- caret::trainControl(
    method = "repeatedcv",
    number = nfold,
    repeats = nfold
  )
  enetFit <- caret::train(
    x = train.x, y = factor(train.y),
    method = "glmnet",
    family = "binomial",
    trControl = fitControl,
    metric = "Accuracy", tuneGrid = grid
  )
  chose_alpha <- enetFit$bestTune[, 1]
  chose_lambda <- enetFit$bestTune[, 2]
  list(chose_alpha = chose_alpha, chose_lambda = chose_lambda)
}



#' Calculate Area Under the Curve (AUC) for Binomial Model
#'
#' This function computes the AUC for a binomial model's predictions on a given dataset.
#' It uses the specified regularization strengths to generate predictions, which are then
#' evaluated against actual outcomes to compute the AUC using the ROCR package. This function
#' is typically used to assess model performance in classification tasks.
#'
#' @param model A model object fitted using a binomial distribution, from which predictions will be generated.
#' @param newx A matrix or data frame of new data on which to make predictions. This should correspond
#'        to the predictors used in fitting the model.
#' @param s A character vector indicating the specific regularization strengths ('lambda.min' or 'lambda.1se')
#'        at which predictions should be evaluated.
#' @param acture.y A vector containing the actual binary outcomes associated with `newx`.
#'
#' @return A numeric value representing the AUC for the model's predictions against actual outcomes.
#'
#' @examples
#' # Assuming 'model', 'newx', and 'actual.y' are predefined:
#' auc_value <- BinomialAUC(model = fitted_model, newx = test_data, s = "lambda.min", acture.y = test_outcomes)
#' print(auc_value)
#' @export
BinomialAUC <- function(model, newx, s, acture.y) {
  prob <- stats::predict(model, newx = newx, s = s, type = "response")
  pred <- ROCR::prediction(prob, acture.y)
  auc <- as.numeric(ROCR::performance(pred, "auc")@y.values)
  return(auc)
}

#' Plot AUC ROC Curves
#'
#' This function generates ROC (Receiver Operating Characteristic) curves to evaluate the performance of a binary classification model.
#' It creates plots for both training and testing data sets, using specified lambda values for model evaluation. The function
#' automatically handles the selection of color palettes and saves the resulting plots in a specified directory.
#'
#' @param train.x A matrix or data frame containing the training predictors.
#' @param train.y A numeric or binary vector indicating the outcomes for the training set.
#' @param test.x A matrix or data frame containing the testing predictors.
#' @param test.y A numeric or binary vector indicating the outcomes for the testing set.
#' @param model A model object used to generate prediction probabilities for the ROC analysis.
#' @param modelname A string representing the name of the model, used in the plot title and file name.
#' @param cols Optional, a vector of colors for the ROC curves. If NULL, default color palettes are applied based on the 'palette' parameter.
#' @param palette Optional, a string specifying the color palette. Default is "jama", which is applied if 'cols' is NULL.
#'
#' @return A ggplot object of the ROC curve plot, which is also saved as a PDF in the specified directory.
#' @examples
#' # Assuming 'train.x', 'train.y', 'test.x', 'test.y', and 'model' are predefined:
#' PlotAUC(
#'   train.x = train_data, train.y = train_outcomes, test.x = test_data, test.y = test_outcomes,
#'   model = fitted_model, modelname = "MyModel"
#' )
#' @export
PlotAUC <- function(train.x, train.y, test.x, test.y, model, modelname, cols = NULL, palette = "jama") {
  if (is.null(cols)) {
    cols <- IOBR::palettes(category = "box", palette = palette, show_message = FALSE, show_col = FALSE)
  } else {
    cols <- cols
  }

  newx <- list(train.x, train.x, test.x, test.x)
  s <- list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y <- list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)
  pref <- args %>% purrr::pmap(CalculatePref, model = model)
  aucs <- args %>%
    purrr::pmap_dbl(BinomialAUC, model = model) %>%
    round(., 2)
  legend.name <- paste(c(
    "train_lambda.min", "train_lambda.1se",
    "test_lambda.min", "test_lambda.1se"
  ), "AUC", aucs, sep = " ")
  names(pref) <- c(
    "train_lambda.min", "train_lambda.1se",
    "test_lambda.min", "test_lambda.1se"
  )
  plotdat <- lapply(pref, function(z) {
    data.frame(x = z@x.values[[1]], y = z@y.values[[1]])
  }) %>% plyr::ldply(., .fun = "rbind", .id = "s")
  plotdat$s <- factor(plotdat$s, levels = names(pref))
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
    ggtitle(str_replace(modelname, "_", " ")) +
    theme(legend.title = element_blank()) +
    theme(
      plot.title = element_text(size = rel(2), hjust = 0.5),
      axis.text.x = element_text(face = "plain", angle = 0, hjust = 1, color = "black"),
      axis.text.y = element_text(face = "plain", angle = 30, hjust = 1, color = "black")
    )

  return(p)
}


#' Calculate Performance Metrics
#'
#' Computes performance metrics such as true positive rate (TPR) and false positive rate (FPR) for model predictions.
#' This function uses the ROCR package to evaluate model effectiveness at different thresholds, helping to
#' assess the discriminative ability of the model under various regularization strengths specified by 's'.
#'
#' @param model A model object used to generate predictions.
#' @param newx A matrix or data frame of new data on which the model will predict.
#' @param s A character vector or single character string indicating the regularization strength at which
#'          the model should be evaluated. Common values are 'lambda.min' and 'lambda.1se' for models fitted with
#'          methods such as glmnet.
#' @param acture.y A vector containing the actual binary outcomes (0 or 1) corresponding to `newx`.
#'
#' @return A performance object from the ROCR package, which includes true positive and false positive rates.
#'
#' @examples
#' # Assuming 'model', 'new_data', and 'actual_outcomes' are predefined:
#' perf_metrics <- CalculatePref(model = fitted_model, newx = new_data, s = "lambda.min", acture.y = actual_outcomes)
#' print(perf_metrics)
#' @export
CalculatePref <- function(model, newx, s, acture.y) {
  prob <- stats::predict(model, newx = newx, s = s, type = "response")
  pred <- ROCR::prediction(prob, acture.y)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  return(perf)
}


#' Split Data into Training and Testing Sets
#'
#' This function divides the dataset into training and testing sets based on a specified proportion. It is designed to work with both binomial and survival type data, ensuring appropriate handling of the dataset
#' according to the specified data type.
#'
#' @param x A matrix or data frame containing the predictors or features.
#' @param y A vector or matrix containing the outcomes associated with 'x'. The format can be either
#'          a simple vector for binomial outcomes or a matrix for survival outcomes.
#' @param train_ratio A numeric value between 0 and 1 indicating the proportion of the data to use for training.
#'                    For example, a 'train_ratio' of 0.7 means 70 percent of the data will be used for training.
#' @param type A character string indicating the type of analysis, accepts "binomial" for binary outcomes or
#'             "survival" for survival analysis, which affects how 'y' data is handled during splitting.
#' @param seed An integer used to set the seed for random sampling, ensuring reproducibility.
#'
#' @return A list containing the training and testing subsets: 'train.x' and 'train.y' for training,
#'         and 'test.x' and 'test.y' for testing. Also returns 'train_sample', which is a vector of indices
#'         that were selected for the training set.
#'
#' @examples
#' # Example for binomial data
#' data_matrix <- matrix(rnorm(200), ncol = 2)
#' outcome_vector <- rbinom(100, 1, 0.5)
#' split_data <- SplitTrainTest(x = data_matrix, y = outcome_vector, train_ratio = 0.7, type = "binomial", seed = 123)
#' print(split_data)
#'
#' # Example for survival data
#' survival_outcomes <- matrix(c(rep(1, 100), rep(0, 100)), ncol = 2)
#' split_data_survival <- SplitTrainTest(x = data_matrix, y = survival_outcomes, train_ratio = 0.7, type = "survival", seed = 123)
#' print(split_data_survival)
#' @export
SplitTrainTest <- function(x, y, train_ratio, type, seed) {
  sizes <- round(nrow(x) * train_ratio)
  set.seed(seed)
  train_sample <- sample(1:nrow(x), size = sizes, replace = F)
  test_sample <- setdiff(1:nrow(x), train_sample)
  train.x <- x[train_sample, ]
  test.x <- x[test_sample, ]
  if (type == "binomial") {
    train.y <- y[train_sample]
    test.y <- y[test_sample]
  }
  if (type == "survival") {
    train.y <- y[train_sample, ]
    test.y <- y[test_sample, ]
  }
  return(list(
    train.x = train.x, train.y = train.y, test.x = test.x, test.y = test.y,
    train_sample = train_sample
  ))
}
