#' Prognostic Model construction
#'
#' @param x input matrix or data.frame where samples are in rows and features in columns; The first column of x is the sample ID of which column is "ID".
#' @param y input matrix or data.frame with three column. Column names are "ID", "time", "status"
#' @param scale A logistic: should the x be scaled, default is TRUE.
#' @param seed default 123456
#' @param train_ratio Value between 0-1; The ratio was used to split the x and y into training and testing data.
#' @param nfold nfold default 10
#' @param plot A logistic, default is TRUE.
#'
#' @return a list contain the results of 2 model (Lasso, Ridge) and the input train data.
#' @export
#'
#' @examples
PrognosticModel <- function(x, y, scale = F, seed = "123456", train_ratio = 0.7, nfold = 10, plot = T){

  x<-as.data.frame(x)
  y<-as.data.frame(y)

  print(message(paste0("\n", ">>> Processing data")))
  processdat <- ProcessingData(x = x, y = y, scale = scale, type = "survival")
  x_scale <- processdat$x_scale
  y <- processdat$y
  x_ID <- processdat$x_ID

  print(message(paste0("\n", ">>> Spliting data into train and test data")))
  train_test <- SplitTrainTest(x = x_scale, y = y, train_ratio = train_ratio,
                               type = "survival", seed = seed)
  train.x = train_test$train.x; train.y <- train_test$train.y
  test.x = train_test$test.x; test.y <- train_test$test.y
  train_sample <- train_test$train_sample
  return.x <- data.frame(ID = x_ID[train_sample], train.x)

  print(message(paste0("\n", ">>> Running ", "LASSO")))
  set.seed(seed)
  lasso_model <- glmnet::cv.glmnet(x = train.x, y = as.matrix(train.y),
                                   family = "cox", alpha = 1, nfolds = nfold)
  lasso_result <- PrognosticResult(model = lasso_model, train.x, train.y, test.x, test.y)
  if (plot){
   p1 <- PlotTimeROC(train.x = train.x, train.y = train.y,
            test.x = test.x, test.y = test.y, model = lasso_model,
            foldername = "5-1_Prognostic_Model",
            modelname = "lasso_model")
   print(p1)
  }

  message(paste0("\n", ">>> Running ", "RIDGE REGRESSION"))

  set.seed(seed)
  ridge_model <- glmnet::cv.glmnet(x = train.x, y = as.matrix(train.y), family = "cox", alpha = 0, nfolds = nfold)
  ridge_result <- PrognosticResult(model = ridge_model, train.x, train.y, test.x, test.y)

  if (plot){
   p2 <- PlotTimeROC(train.x = train.x, train.y = train.y,
                test.x = test.x, test.y = test.y, model = ridge_model,
                foldername = "5-1_Prognostic_Model",
                modelname = "ridge_model")
   print(p2)
  }
  return(list(lasso_result = lasso_result, ridge_result = ridge_result,
              train.x = return.x))
  message(paste0("\n", ">>> Done !"))
}

#' Prognostic Result
#'
#' @param model
#' @param train.x
#' @param train.y
#' @param test.x
#' @param test.y
#'
#' @return
#' @export
#'
#' @examples
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


#' Prognostic AUC
#'
#' @param model
#' @param newx
#' @param s
#' @param acture.y
#'
#' @return
#' @export
#'
#' @examples
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



#' Calculate Time ROC
#'
#' @param model
#' @param newx
#' @param s
#' @param acture.y
#' @param foldername
#' @param modelname
#'
#' @return
#' @export
#'
#' @examples
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

#' Plot Time ROC
#'
#' @param train.x
#' @param train.y
#' @param test.x
#' @param test.y
#' @param model
#' @param foldername
#' @param modelname
#'
#' @return
#' @export
#'
#' @examples
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
  names(roclist) <- c("train_lambda.min", "train_lambda.1se",
                   "test_lambda.min", "test_lambda.1se")

  plotdat <- lapply(roclist, function(z){
    data.frame(x = z$FP[, 2], y = z$TP[, 2])
  }) %>% plyr::ldply(., .fun = "rbind", .id = "s")
  plotdat$s <- factor(plotdat$s, levels = names(roclist))

  p <- ggplot2::ggplot(plotdat, aes(x = x, y = y)) +
    geom_path(aes(color= s)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("False positive rate") + ylab("True positive rate") +
    theme_bw() + scale_color_manual(values = mycols,
                                    labels = legend.name) +
    ggtitle(paste0(str_replace(modelname, "_", " "), "\nROC at time quantile 0.9")) +
    theme(legend.title = element_blank()) +
    theme(plot.title=element_text(size=rel(1.5),hjust=0.5),
          axis.text.x= element_text(face="plain",angle=0,hjust = 1,color="black"),
          axis.text.y= element_text(face="plain",angle=30,hjust = 1,color="black"))
  if (!dir.exists(foldername)){
    dir.create(foldername)}
  ggplot2::ggsave(paste0(foldername, "/", modelname, "_ROC.pdf"), plot = p, width = 6, height = 4)
  return(p)
}

