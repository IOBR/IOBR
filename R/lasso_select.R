

#' Constructing predictive or prognostic modle
#'
#' @param x input matrix.Rownames should be features like gene symbols or cgi, colnames be samples
#' @param y response variable. can be binary and survival type.
#' @param type "binary" or "survival"
#' @param nfold number of nfold for cross-validation - default is 10
#' @param lambda "lambda.min" or "lambda.1se"
#'
#' @return
#' @export
#'
#' @examples
#'
#' data("crc_clin")
#' data("tcga_crc_exp")
#' dat <- t(tcga_crc_exp)
#' dat <- data.frame(patient = str_sub(rownames(dat), 1, 12), dat)
#' dat <- merge(crc_clin, dat, by = "patient")
#' dat <- dat[dat$OS_time > 0, ]
#' x <- t(dat[, -c(1:3)])
#' mad <- apply(x, 1, mad)
#' x <- x[mad > 0.5, ]
#' y <- Surv(dat$OS_time, dat$OS)
#' pd1 <- as.numeric(dat[, "PDCD1"])
#' group <- ifelse(pd1 > mean(pd1), 1, 0)
#' sur_gene <- lasso_select(x = x, y = y, type = "survival", nfold = 10, lambda = "lambda.min")
#' pd1_gene <- lasso_select(x = x, y = group, type = "binary", nfold = 10, lambda = "lambda.min")
lasso_select <- function(x, y, type =c("binary", "survival"), nfold = 10,
                         lambda = c("lambda.min", "lambda.1se")){
  type = match.arg(type)
  lambda = match.arg(lambda)
  x <- t(x)
  if (type == "binary"){
    cvfit = cv.glmnet(x, y,
                      nfold,
                      type.measure = "class")
  }else{
    cvfit = cv.glmnet(x, y, nfold, family = "cox")
  }
  myCoefs <- coef(cvfit, s = lambda);
  lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
  lasso_fea <- lasso_fea[-1]
  feature <- lasso_fea
  return(feature)
}
