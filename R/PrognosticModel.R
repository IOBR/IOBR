#' Title
#'
#' @param x input matrix or data.frame where samples are in rows and features in columns; The first column of x is the sample ID of which column is "ID".
#' @param y input matrix or data.frame with three column. Column names are "ID", "time", "status"
#' @param scale A logistic: should the x be scaled, default is TRUE.
#' @param seed default 123456
#' @param train_ratio Value between 0-1; The ratio was used to split the x and y into training and testing data.
#' @param nfold nfold default 10
#' @param plot
#'
#' @return a list contain the results of 2 model (Lasso, Ridge) and the input train data.
#' @export
#' @import glmnet
#' @import timeROC
#' @import stats
#' @import purrr
#' @import dplyr
#'
#'
#' @examples
#' load("./data/crc_response.rda")
#' load("./data/tme_combine.rda")
#' x <- tme_combine
#' y <- na.omit(response[, c("ID", "OS.time", "OS")])
#' colnames(y) <- c("ID", "time", "status")
#' PrognosticResult <- PrognosticModel(x = x, y = y, scale = T, seed = "123456", train_ratio = 0.7)
PrognosticModel <- function(x, y, scale = T, seed = "123456", train_ratio = 0.7, nfold = 10, plot = T){
  message(paste0("\n", ">>> Processing data"))
  processdat <- ProcessingData(x = x, y = y, scale = scale, type = "survival")
  x_scale <- processdat$x_scale
  y <- processdat$y
  x_ID <- processdat$x_ID

  message(paste0("\n", ">>> Spliting data into train and test data"))
  train_test <- SplitTrainTest(x = x_scale, y = y, train_ratio = train_ratio, type = "survival")
  train.x = train_test$train.x; train.y <- train_test$train.y
  test.x = train_test$test.x; test.y <- train_test$test.y
  train_sample <- train_test$train_sample
  return.x <- data.frame(ID = x_ID[train_sample], train.x)

  message(paste0("\n", ">>> Running ", "LASSO"))
  set.seed(seed)
  lasso_model <- glmnet::cv.glmnet(x = train.x, y = as.matrix(train.y), family = "cox", alpha = 1, nfolds = nfold)
  lasso_result <- PrognosticResult(model = lasso_model, train.x, train.y, test.x, test.y)
  if (plot){
    PlotTimeROC(train.x = train.x, train.y = train.y,
            test.x = test.x, test.y = test.y, model = lasso_model,
            foldername = "5-1_Prognostic_Model",
            modelname = "lasso_model")
  }

  message(paste0("\n", ">>> Running ", "RIDGE REGRESSION"))

  set.seed(seed)
  ridge_model <- glmnet::cv.glmnet(x = train.x, y = as.matrix(train.y), family = "cox", alpha = 0, nfolds = nfold)
  ridge_result <- PrognosticResult(model = ridge_model, train.x, train.y, test.x, test.y)

  if (plot){
    PlotTimeROC(train.x = train.x, train.y = train.y,
                test.x = test.x, test.y = test.y, model = ridge_model,
                foldername = "5-1_Prognostic_Model",
                modelname = "ridge_model")
  }

  if (plot){
      message(paste0("\n", ">>> Generate figures are in the current working directory"))
      models <- list(lasso_model, ridge_model)
      foldername <- list("5-1_Prognostic_Model")
      resultnames <- list("lasso_model", "ridge_model")
      agrs <- list(model = models, foldername = foldername, resultname = resultnames) %>%
        purrr::pmap(PlotModel)
    }

  return(list(lasso_result = lasso_result, ridge_result = ridge_result,
              train.x = return.x))
  message(paste0("\n", ">>> Done !"))
}

PrognosticResult <- function(model, train.x, train.y, test.x, test.y){
  coefs <- cbind(stats::coef(model, s = "lambda.min"), stats::coef(model, s = "lambda.1se"))
  coefs <- data.frame(feature = rownames(coefs), lambda.min =coefs[, 1], lambda.1se = coefs[, 2])
  newx = list(train.x, train.x, test.x, test.x)
  s = list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y = list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)
  auc <-  args %>% purrr::pmap(PrognosticAUC, model = model) %>% do.call(rbind, .)
  rownames(auc) <- c("Train_lambda.min", "Train_lambda.1se", "Test_lambda.min", "Test_lambda.1se")
  resultreturn <- list(model = model, coefs = coefs,
                       AUC = auc)
}
PrognosticAUC <- function(model, newx, s, acture.y){
  riskscore <- stats::predict(model, newx = newx, s = s)
  timerocDat <- data.frame(risk = riskscore[, 1], acture.y)
  with(timerocDat,
       ROC <<- timeROC::timeROC(T = time, delta = status,
                       marker = risk, cause = 1,
                       weighting = "marginal",
                       time = quantile(time,probs = c(0.3, 0.9)),
                       ROC = TRUE,
                       iid = TRUE))
  AUC <- data.frame(probs.3 = ROC$AUC[1], probs.9 = ROC$AUC[2])
  return(AUC)
}

CalculateTimeROC <- function(model, newx, s, acture.y, foldername, modelname){
  riskscore <- stats::predict(model, newx = newx, s = s)
  timerocDat <- data.frame(risk = riskscore[, 1], acture.y)
  with(timerocDat,
       ROC <<- timeROC::timeROC(T = time, delta = status,
                                marker = risk, cause = 1,
                                weighting = "marginal",
                                time = quantile(time,probs =0.9),
                                ROC = TRUE,
                                iid = TRUE))
  ROC
}

PlotTimeROC <- function(train.x, train.y, test.x, test.y, model, foldername, modelname){
  mycols <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
              "#F39B7FFF", "#8491B4FF", "#91D1C2FF")
  newx = list(train.x, train.x, test.x, test.x)
  s = list("lambda.min", "lambda.1se", "lambda.min", "lambda.1se")
  acture.y = list(train.y, train.y, test.y, test.y)
  args <- list(newx, s, acture.y)
  auc <-  args %>% purrr::pmap(PrognosticAUC, model = model) %>% do.call(rbind, .)
  rownames(auc) <- c("Train_lambda.min", "Train_lambda.1se", "Test_lambda.min", "Test_lambda.1se")
  roclist <- args %>% purrr::pmap(CalculateTimeROC, model = model)
  aucs <- round(auc$probs.9, 2)
  legend.name <- paste(c("train_lambda.min", "train_lambda.1se",
                         "test_lambda.min", "test_lambda.1se"), "AUC", aucs, sep=" ")
  if (!dir.exists(foldername)){
    dir.create(foldername)}
  pdf(paste0(foldername, "/", modelname, "_ROC.pdf"), width = 6, height = 6)
  plot(roclist[[1]], col = mycols[1], add = FALSE,
       time = as.numeric(roclist[[1]]$times)[2], title = FALSE)
  for (i in 2:length(roclist)){
    plot(roclist[[i]], col = mycols[i], add = TRUE,
         time = as.numeric(roclist[[i]]$times)[2])
  }
  title("ROC at time quantile 0.9")
  legend("bottomright", legend = legend.name,
         col = mycols[1:length(roclist)], lwd = 2)
  dev.off()
}

