#' Title
#'
#' @param x data.frame contains sample ID and features; The first column of x is the sample ID.
#' @param y data.frame whose sample ID in the first column and the outcome of each sample in the second sample. Outcome value can be numeric or factor vector.
#' @param seed default 123456
#' @param scale A logistic: should the x be scaled, default is TRUE.
#' @param train_ratio Value between 0-1, eg: 0.7; The ratio is used to split the x and y into training and testing data.
#' @param nfold default 10
#' @param plot
#'
#' @return a list contain the results of 3 model (Lasso, Ridge and ElasticNet regression model) and the input train data.
#'
#' @export
#'
#' @import glmnet
#' @import caret
#' @import pROC
#' @import dplyr
#' @import magrittr
#'
#' @examples
#' load("./data/crc_response.rda")
#' load("./data/tme_combine.rda")
#' binomial_result <- BinomialModel(x = tme_combine, y = na.omit(response[, c("ID", "Stage")]),
#' seed = "123456", scale = TRUE, train_ratio = 0.7, nfold = 10)
BinomialModel <- function(x, y,seed = "123456", scale = TRUE,
                          train_ratio = 0.7, nfold = 10, plot = F){
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
  lasso_model <- glmnet::cv.glmnet(x = train.x, y = train.y, family = "binomial",
                     type.measure = "class", alpha = 1,  nfolds = nfold)
  lasso_result <- RegressionResult(train.x = train.x, train.y = train.y,
                                   test.x = test.x, test.y = test.y, model = lasso_model)
  if (plot){
    PlotAUC(train.x = train.x, train.y = train.y,
            test.x = test.x, test.y = test.y, model = lasso_model,
            foldername = "5-1_Binomial_Model",
            modelname = "lasso_model")
  }

  message(paste0("\n", ">>> Running ", "RIDGE REGRESSION"))

  set.seed(seed)
  ridge_model <- glmnet::cv.glmnet(x = train.x, y = train.y, family = "binomial",
                           type.measure = "class", alpha = 0, nfolds = nfold)
  ridge_result <- RegressionResult(train.x = train.x, train.y = train.y,
                                   test.x = test.x, test.y = test.y, model = ridge_model)
  if (plot){
    PlotAUC(train.x = train.x, train.y = train.y,
            test.x = test.x, test.y = test.y, model = ridge_model,
            foldername = "5-1_Binomial_Model",
            modelname = "ridge_model")
  }
  message(paste0("\n", ">>> Running ", "Elastic Network."))


  lambdamax <- max(lasso_model$lambda.min, lasso_model$lambda.1se,
                   ridge_model$lambda.min, ridge_model$lambda.1se)
  lambdamax <- round(lambdamax, 1) + 0.1
  if (lambdamax > 1){lambdamax <- 1}
  message(paste0("\n", ">>> Choosing optimal alpha and lambda.", "\n This step may take some time."))

  set.seed(seed)
  AlphaLambda <- Enet(train.x, train.y, lambdamax, nfold = nfold)
  set.seed(seed)
  enet_model <- glmnet::glmnet(x = train.x, y = train.y, family = "binomial",
                          type.measure = "class",
                          alpha = AlphaLambda$chose_alpha, lambda = AlphaLambda$chose_lambda)
  enet_result <- RegressionResult(train.x = train.x, train.y = train.y,
                                  test.x = test.x, test.y = test.y, model = enet_model)
  if (plot){
    PlotAUC(train.x = train.x, train.y = train.y,
            test.x = test.x, test.y = test.y, model = lasso_model,
            foldername = "5-1_Binomial_Model",
            modelname = "enet_model")
  }

  enet_result <- c(enet_result, AlphaLambda)
  if (plot){
    message(paste0("\n", ">>> Generate figures are in the current working directory"))
    models <- list(lasso_model, ridge_model, enet_model)
    foldername <- list("5-1_Binomial_Model")
    resultnames <- list("lasso_model", "ridge_model", "enet_model")
    agrs <- list(model = models, foldername = foldername, resultname = resultnames) %>%
      purrr::pmap(PlotModel)
  }

  return(list(lasso_result = lasso_result, ridge_result = ridge_result,
              enet_result = enet_result, train.x = return.x))

  message(paste0("\n", ">>> Done !"))
}


#######################################################
ProcessingData <- function(x, y, scale = scale, type = "binomial"){
  colnames(x)[1] <- "ID"
  colnames(y)[1] <- "ID"
  x$ID <- as.character(x$ID)
  y$ID <- as.character(y$ID)
  if (type == "survival"){
    colnames(y) <- c("ID", "time", "status")
    y <- dplyr::filter(y, time > 0)
  }
  samples <- intersect(x$ID, y$ID)

  if (length(samples) == 0){
    stop("No same sample ID been found between input matrix x and y")
  }
  x <- x[match(samples, x$ID), ]
  y <- y[match(samples, y$ID), ]
  if (type == "binomial"){
    colnames(y) <- c("ID", "Group")
    y <- dplyr::pull(y, Group)
    if (!is.factor(y)){
      message(paste0("\n", ">>> Outcome is not a factor, transform it into factor vector."))
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
  if (length(ValueNA) > 0){x_scale <- x_scale[, -ValueNA]}


  return(list(x_scale = x_scale, y = y, x_ID = x_ID))
}


RegressionResult <- function(train.x, train.y, test.x, test.y, model){
  coefs <- cbind(coef(model, s = "lambda.min"), coef(model, s = "lambda.1se"))
  coefs <- data.frame(feature = rownames(coefs), lambda.min =coefs[, 1], lambda.1se = coefs[, 2])
  newx = list(train.x, train.x, test.x, test.x)
  s = list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y = list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)
  AUC <-  args %>% purrr::pmap_dbl(BinomialAUC, model = model) %>%
    matrix(., ncol = 2, byrow = T,
           dimnames = list(c("train", "test"), c("lambda.min", "lambda.1se")))
  resultreturn <- list(model = model, coefs = coefs,
                       AUC = AUC)
}
Enet <- function(train.x, train.y, lambdamax, nfold = nfold){
  grid <- expand.grid(.alpha = seq(0, 1, by = .2), .lambda = seq(0, lambdamax, length.out = 10))
  fitControl <- caret::trainControl(method = "repeatedcv",
                             number = nfold,
                             repeats = nfold)
  enetFit <- caret::train(x = train.x, y = factor(train.y),
                     method = "glmnet",
                     family = "binomial",
                     trControl = fitControl,
                     metric = "Accuracy", tuneGrid = grid)
  chose_alpha <- enetFit$bestTune[, 1]
  chose_lambda <- enetFit$bestTune[, 2]
  list(chose_alpha = chose_alpha, chose_lambda = chose_lambda)
  }

BinomialAUC <- function(model, newx, s, acture.y){
  prob <- stats::predict(model, newx = newx, s = s, type = "response")
  pred <- ROCR::prediction(prob, acture.y)
  auc <- as.numeric(ROCR::performance(pred, "auc")@y.values)
  return(auc)
}
PlotAUC <- function(train.x, train.y, test.x, test.y, model, foldername, modelname){
  mycols <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
              "#F39B7FFF", "#8491B4FF", "#91D1C2FF")
  newx = list(train.x, train.x, test.x, test.x)
  s = list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y = list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)
  pref <- args %>% purrr::pmap(CalculatePref, model = model)
  aucs <- args %>% purrr::pmap_dbl(BinomialAUC, model = model) %>% round(., 2)
  legend.name <- paste(c("train_lambda.min", "train_lambda.1se",
                         "test_lambda.min", "test_lambda.1se"), "AUC", aucs,sep=" ")
  if (!dir.exists(foldername)){
    dir.create(foldername)}
  pdf(paste0(foldername, "/", modelname, "_ROC.pdf"), width = 6, height = 6)
  plot(pref[[1]], main = "ROC", col = mycols[1])
  for (i in 2:length(pref)){
    plot(pref[[i]], main = "ROC", col = mycols[i], add = T)
  }
  legend("bottomright", legend = legend.name,
         col = mycols[1:4], lwd = 2)
  dev.off()
}
CalculatePref<- function(model, newx, s, acture.y){
  prob <- stats::predict(model, newx = newx, s = s, type = "response")
  pred <- ROCR::prediction(prob, acture.y)
  perf <- ROCR::performance(pred,"tpr","fpr")
  return(perf)
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

PlotModel <- function(model, foldername, resultname){
  if (!dir.exists(foldername)){
    dir.create(foldername)}
    pdf(paste0(foldername, "/", resultname, ".pdf"), width = 6, height = 4)
    plot(model)
    dev.off()
  }


