#' Extract features of Predictive or Prognostic Model Using LASSO Regression
#'
#' This function applies LASSO (Least Absolute Shrinkage and Selection Operator) regression to construct predictive or
#' prognostic models. It is designed to handle both binary and survival type response variables, using cross-validation
#' to optimize the model selection.
#'
#' @param x An input matrix with features like gene symbols or CGI as row names and samples as column names.
#' @param y A response variable vector; can be binary (0 or 1) or survival data (survival time and event status).
#' @param type A character string specifying the model type: "binary" for binary response models or "survival" for survival analysis.
#'        The default is set to "binary".
#' @param nfold The number of folds for cross-validation, with a default of 10.
#' @param lambda The regularization parameter selection method: "lambda.min" for the lambda that gives minimum mean
#'        cross-validated error or "lambda.1se" for the lambda that is one standard error away from the minimum.
#'        The default is set to "lambda.min".
#'
#' @return A vector of selected feature names from the LASSO model that are non-zero in the optimal model.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("gene_expression", package = "examplePackage")
#' survival_data <- survival::Surv(time = c(2, 4, 3), event = c(1, 0, 1))
#' # Example for a binary response variable
#' binary_outcome <- c(0, 1, 1)
#' model_features_binary <- lasso_select(x = gene_expression, y = binary_outcome, type = "binary")
#' # Example for a survival response variable
#' model_features_survival <- lasso_select(x = gene_expression, y = survival_data, type = "survival")
lasso_select <- function(x, y, type = c("binary", "survival"), nfold = 10,
                         lambda = c("lambda.min", "lambda.1se")) {
  type <- match.arg(type)
  lambda <- match.arg(lambda)
  x <- t(x)
  if (type == "binary") {
    cvfit <- cv.glmnet(x, y,
      nfold,
      alpha = 1,
      type.measure = "class"
    )
  } else {
    cvfit <- cv.glmnet(x, y, nfold, alpha = 1, family = "cox")
  }
  myCoefs <- coef(cvfit, s = lambda)
  lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0)]
  lasso_fea <- lasso_fea[-1]
  feature <- lasso_fea
  return(feature)
}
