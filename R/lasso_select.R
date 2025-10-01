#' Feature Selection for Predictive or Prognostic Models Using LASSO Regression
#'
#' Applies LASSO (Least Absolute Shrinkage and Selection Operator) regression to construct predictive or prognostic models. Supports both binary and survival response variables, utilizing cross-validation for optimal model selection.
#'
#' @param x A numeric matrix. Features (e.g., gene symbols or CGI) as row names and samples as column names.
#' @param y A response variable vector. Can be binary (0/1) or survival data (e.g., survival time and event status).
#' @param type Character. Model type: "binary" for binary response or "survival" for survival analysis. Default is "binary".
#' @param nfold Integer. Number of folds for cross-validation. Default is 10.
#' @param lambda Character. Regularization parameter selection: "lambda.min" (minimum mean cross-validated error) or "lambda.1se" (one standard error from minimum). Default is "lambda.min".
#'
#' @return Character vector of selected feature names with non-zero coefficients in the optimal LASSO model.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("gene_expression", package = "examplePackage")
#' survival_data <- survival::Surv(time = c(2, 4, 3), event = c(1, 0, 1))
#' # Binary response example
#' binary_outcome <- c(0, 1, 1)
#' lasso_select(x = gene_expression, y = binary_outcome, type = "binary")
#' # Survival response example
#' lasso_select(x = gene_expression, y = survival_data, type = "survival")
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
