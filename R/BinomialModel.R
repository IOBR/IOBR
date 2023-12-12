#' Binomial Model construction
#'
#' @param x data.frame contains sample ID and features; The first column of x is the sample ID.
#' @param y data.frame whose sample ID in the first column and the outcome of each sample in the second sample. Outcome value can be numeric or factor vector.
#' @param seed default 123456
#' @param scale A logistic: should the x be scaled, default is TRUE.
#' @param train_ratio Value between 0-1, eg: 0.7; The ratio is used to split the x and y into training and testing data.
#' @param nfold default 10
#' @param plot A logistic, default is TRUE.
#'
#' @return a list contain the results of 2 model (Lasso, Ridge) and the input train data.
#'
#' @export
#' @examples
#' data("imvigor210_sig", package = "IOBR")
#' data("imvigor210_pdata",package = "IOBR")
#' pdata_group <- imvigor210_pdata[!imvigor210_pdata$BOR_binary=="NA",c("ID","BOR_binary")]
#' pdata_group$BOR_binary <- ifelse(pdata_group$BOR_binary == "R", 1, 0)
#' BinomialModel(x = imvigor210_sig, y = pdata_group, seed = 123456, scale = TRUE, train_ratio = 0.7, nfold = 10, plot = T)
BinomialModel <- function(x, y,seed = 123456, scale = TRUE, train_ratio = 0.7, nfold = 10, plot = T){

  x<-as.data.frame(x)
  y<-as.data.frame(y)
  print(message(paste0("\n", ">>> Processing data")))
  processdat <- ProcessingData(x = x, y = y, scale = scale, type = "binomial")
  x_scale <- processdat$x_scale
  y <- processdat$y
  x_ID <-processdat$x_ID
  print(message(paste0("\n", ">>> Spliting data into train and test data")))
  train_test <- SplitTrainTest(x = x_scale, y = y, train_ratio = train_ratio, type = "binomial",
                               seed = seed)
  train.x = train_test$train.x; train.y <- train_test$train.y
  test.x = train_test$test.x; test.y <- train_test$test.y
  train_sample <- train_test$train_sample
  return.x <- data.frame(ID = x_ID[train_sample], train.x)

  print(message(paste0("\n", ">>> Running ", "LASSO")))
  set.seed(seed)
  lasso_model <- glmnet::cv.glmnet(x = train.x, y = train.y, family = "binomial",
                     type.measure = "class", alpha = 1,  nfolds = nfold)
  lasso_result <- RegressionResult(train.x = train.x, train.y = train.y,
                                   test.x = test.x, test.y = test.y, model = lasso_model)
  if (plot){
    p1 <- PlotAUC(train.x = train.x, train.y = train.y,
                 test.x = test.x, test.y = test.y, model = lasso_model,
                 foldername = "5-1_Binomial_Model",
                 modelname = "lasso_model")
    print(p1)
  }

  print(message(paste0("\n", ">>> Running ", "RIDGE REGRESSION")))

  set.seed(seed)
  ridge_model <- glmnet::cv.glmnet(x = train.x, y = train.y, family = "binomial",
                           type.measure = "class", alpha = 0, nfolds = nfold)
  ridge_result <- RegressionResult(train.x = train.x, train.y = train.y,
                                   test.x = test.x, test.y = test.y, model = ridge_model)
  if (plot){
    p2 <- PlotAUC(train.x = train.x, train.y = train.y,
            test.x = test.x, test.y = test.y, model = ridge_model,
            foldername = "5-1_Binomial_Model",
            modelname = "ridge_model")
    print(p2)
  }
  print(message(paste0("\n", ">>> Running ", "Elastic Network.")))

  return(list(lasso_result = lasso_result, ridge_result = ridge_result,
              train.x = return.x))

  message(paste0("\n", ">>> Done !"))
}


#######################################################
#' Processing Data
#'
#' @param x
#' @param y
#' @param scale
#' @param type
#'
#' @return
#' @export
#'
#' @examples
ProcessingData <- function(x, y, scale, type = "binomial"){
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


#' Regression Result
#'
#' @param train.x
#' @param train.y
#' @param test.x
#' @param test.y
#' @param model
#'
#' @return
#' @export
#'
#' @examples
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

#' Enet
#'
#' @param train.x
#' @param train.y
#' @param lambdamax
#' @param nfold
#'
#' @return
#' @export
#'
#' @examples
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



#' BinomialAUC
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
BinomialAUC <- function(model, newx, s, acture.y){
  prob <- stats::predict(model, newx = newx, s = s, type = "response")
  pred <- ROCR::prediction(prob, acture.y)
  auc <- as.numeric(ROCR::performance(pred, "auc")@y.values)
  return(auc)
}

#' Plot AUC
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
  names(pref) <- c("train_lambda.min", "train_lambda.1se",
                   "test_lambda.min", "test_lambda.1se")
  plotdat <- lapply(pref, function(z){
    data.frame(x = z@x.values[[1]], y = z@y.values[[1]])
  }) %>% plyr::ldply(., .fun = "rbind", .id = "s")
  plotdat$s <- factor(plotdat$s, levels = names(pref))
  p <- ggplot2::ggplot(plotdat, aes(x = x, y = y)) +
    geom_path(aes(color= s)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("False positive rate") + ylab("True positive rate") +
    theme_bw() + scale_color_manual(values = mycols,
                                    labels = legend.name) +
    ggtitle(str_replace(modelname, "_", " ")) +
    theme(legend.title = element_blank()) +
    theme(plot.title=element_text(size=rel(2),hjust=0.5),
            axis.text.x= element_text(face="plain",angle=0,hjust = 1,color="black"),
            axis.text.y= element_text(face="plain",angle=30,hjust = 1,color="black"))

  if (!dir.exists(foldername)){
    dir.create(foldername)}
  ggplot2::ggsave(paste0(foldername, "/", modelname, "_ROC.pdf"), plot = p, width = 6, height = 4)
  return(p)
}


#' Calculate Pref
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
CalculatePref<- function(model, newx, s, acture.y){
  prob <- stats::predict(model, newx = newx, s = s, type = "response")
  pred <- ROCR::prediction(prob, acture.y)
  perf <- ROCR::performance(pred,"tpr","fpr")
  return(perf)
}

#' Split Train and Test data
#'
#' @param x
#' @param y
#' @param train_ratio
#' @param type
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
SplitTrainTest <- function(x, y, train_ratio, type, seed){
  sizes <- round(nrow(x) * train_ratio)
  set.seed(seed)
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


