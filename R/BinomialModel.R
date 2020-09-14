#' Title
#'
#' @param x input matrix or data.frame where samples are in rows and features in columns; The first column of x is the sample ID of which column is "ID".
#' @param y input matrix or data.frame whose sample ID in the first column and the outcome of each sample in the second sample. Outcome value can be numeric or factor vector.
#' @param seed default 123456
#' @param scale A logistic: should the x be scaled, default is TRUE.
#' @param train_ratio Value between 0 -1; The ratio was used to split the x and y into training and testing data.
#' @param nfold
#'
#' @return a list contain the results of 3 model (Lasso, Ridge and ElasticNet regression model) and the input train data.
#'
#' @export
#'
#' @import glmnet
#' @import caret
#' @import pROC
#'
#' @examples
#' load("./data/crc_response.rda")
#' load("./data/tme_combine.rda")
#' binomial_result <- BinomialModel(x = tme_combine, y = na.omit(response[, c("ID", "Stage")]),
#' seed = "123456", scale = TRUE, train_ratio = 0.7, nfold = 10)
BinomialModel <- function(x, y,seed = "123456", scale = TRUE, train_ratio = 0.7, nfold = 10){
  message(paste0("\n", ">>> Processing data"))
  processdat <- ProcessingData(x = x, y = y, scale = scale, type = "binomial")
  x_scale <- processdat$x_scale
  y <- processdat$y
  x_ID <-processdat$x_ID
  message(paste0("\n", ">>> Spliting data into train and test data"))
  train_test <- SplitTrainTest(x = x_scale, y = y, train_ratio = train_ratio, type = "binomial")
  train.x = train_test$train.x; train.y <- train_test$train.y
  test.x = train_test$test.x; test.y <- train_test$test.y
  train_sample <- train_test$train_sample
  return.x <- data.frame(ID = x_ID[train_sample], train.x)

  message(paste0("\n", ">>> Running ", "LASSO"))
  set.seed(seed)
  lasso_model <- cv.glmnet(x = train.x, y = train.y, family = "binomial",
                     type.measure = "class", alpha = 1,  nfolds = nfold)
  lasso_result <- RegressionResult(train.x = train.x, train.y = train.y,
                                   test.x = test.x, test.y = test.y, model = lasso_model)

  message(paste0("\n", ">>> Running ", "RIDGE REGRESSION"))

  set.seed(seed)
  ridge_model <- cv.glmnet(x = train.x, y = train.y, family = "binomial",
                           type.measure = "class", alpha = 0, nfolds = nfold)
  ridge_result <- RegressionResult(train.x = train.x, train.y = train.y,
                                   test.x = test.x, test.y = test.y, model = ridge_model)

  message(paste0("\n", ">>> Running ", "Elastic Network."))


  lambdamax <- max(lasso_model$lambda.min, lasso_model$lambda.1se,
                   ridge_model$lambda.min, ridge_model$lambda.1se)
  lambdamax <- round(lambdamax, 1) + 0.1
  if (lambdamax > 1){lambdamax <- 1}
  message(paste0("\n", ">>> Choosing optimal alpha and lambda.", "\n This step may take some time."))

  set.seed(seed)
  AlphaLambda <- Enet(train.x, train.y, lambdamax, nfold = nfold)
  set.seed(seed)
  enet_model <- glmnet(x = train.x, y = train.y, family = "binomial",
                          type.measure = "class",
                          alpha = AlphaLambda$chose_alpha, lambda = AlphaLambda$chose_lambda)
  enet_result <- RegressionResult(train.x = train.x, train.y = train.y,
                                  test.x = test.x, test.y = test.y, model = enet_model)
  enet_result <- c(enet_result, AlphaLambda)
  return(list(lasso_result = lasso_result, ridge_result = ridge_result,
              enet_result = enet_result, train.x = return.x))

  message(paste0("\n", ">>> Done !"))
}


#######################################################
ProcessingData <- function(x, y, scale = scale, type = "binomial"){
  x$ID <- as.character(x$ID)
  y$ID <- as.character(y$ID)
  if (type == "survival"){
    colnames(y) <- c("ID", "time", "status")
    y <- dplyr::filter(y, time > 0)
  }
  samples <- intersect(x$ID, y$ID)

  if (length(samples) == 0){
    stop("No same ID been found between input matrix x and y")
  }
  x <- x[match(samples, x$ID), ]
  y <- y[match(samples, y$ID), ]
  if (type == "binomial"){
    colnames(y) <- c("ID", "Group")
    y <- pull(y, Group)
    if (is.character(y)){
      message(paste0("\n", ">>> Outcome value in y is a character, transform it into factor vector."))
      y <- y %>% as.factor()
    }
  }
  if (type == "survival"){
    y <- y[, c("time", "status")]
  }
  if (scale){
    x_scale <- scale(x[, -1], center = TRUE, scale = TRUE)
  }else{x_scale <- x[, -1]}
  x_ID <- x[, "ID"]
  ValueNA <- which(as.numeric(apply(x_scale, 2, function(z)sum(is.na(z)))) != 0)
  x_scale <- x_scale[, -ValueNA]

  return(list(x_scale = x_scale, y = y, x_ID = x_ID))
}


RegressionResult <- function(train.x, train.y, test.x, test.y, model){
  coefs <- cbind(coef(model, s = "lambda.min"), coef(model, s = "lambda.1se"))
  coefs <- data.frame(feature = rownames(coefs), lambda.min =coefs[, 1], lambda.1se = coefs[, 2])
  newx = list(train.x, train.x, test.x, test.x)
  s = list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y = list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)
  AUC <-  args %>% pmap_dbl(BinomialAUC, model = model) %>%
    matrix(., ncol = 2, byrow = T, dimnames = list(c("train", "test"), c("lambda.min", "lambda.1se")))
  resultreturn <- list(model = model, coefs = coefs,
                       AUC = AUC)
}
Enet <- function(train.x, train.y, lambdamax, nfold = nfold){
  grid <- expand.grid(.alpha = seq(0, 1, by = .2), .lambda = seq(0, lambdamax, length.out = 10))
  fitControl <- trainControl(method = "repeatedcv",
                             number = nfold,
                             repeats = nfold)
    enetFit <- train(x = train.x, y = factor(train.y),
                     method = "glmnet",
                     family = "binomial",
                     trControl = fitControl,
                     metric = "Accuracy", tuneGrid = grid)
  chose_alpha <- enetFit$bestTune[, 1]
  chose_lambda <- enetFit$bestTune[, 2]
  list(chose_alpha = chose_alpha, chose_lambda = chose_lambda)
  }

BinomialAUC <- function(model, newx, s, acture.y){
  pred <- predict(model, newx = newx, s = s, type = "class")
  auc <- as.numeric(pROC::roc(predictor = as.numeric(pred[, 1]), response = acture.y)$auc)
  return(auc)
}

SplitTrainTest <- function(x, y, train_ratio = 0.7, type){
  sizes <- round(nrow(x) * train_ratio)
  train_sample <- sample(1:nrow(x), size = sizes, replace = F)
  test_sample <- setdiff(1:nrow(x), train_sample)
  train.x <- x[train_sample, ];
  test.x <- x[test_sample, ];
  if (type == "binomial"){
    train.y <- y[train_sample]
    test.y <- y[test_sample]
  }
  if (type == "survival"){
    train.y <- y[train_sample, ]
    test.y <- y[test_sample, ]
  }
  return(list(train.x = train.x, train.y = train.y, test.x = test.x, test.y = test.y,
              train_sample = train_sample))
}
