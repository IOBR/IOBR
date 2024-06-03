
#' Feature selection
#'
#' This function selects features based on correlation or differential expression methods. 
#' It integrates tools from the limma package for differential expression and basic statistical tests for correlation analysis.
#'
#' @param x input matrix.Rownames should be features like gene symbols or cgi, colnames be samples
#' @param y response variable. Data type can be quantitative, binary and survival. Survival type can be generated through ?survival::Surv
#' @param family For method="cor", useser can choose "spearman" or "pearson" .
#' @param method Binary for method = "dif", quantitative response value for "dif" and "cor".
#' @param cutoff Numeric. Estimate and log2FC cutoff value for correlation analysis and limma dif analysis.
#' @param padjcut Numeric. Adjust P value cutoff.
#'
#' @import tidyverse
#' @import glmnet
#' @import limma
#'
#' @return Returns a vector of selected feature names based on the specified criteria.
#' @export
#'
#' @examples
#'
#'data("imvigor210_eset",package = "IOBR")
#'data("imvigor210_pdata", package = "IOBR")
#'
#'mad <- apply(imvigor210_eset, 1, mad)
#'imvigor210_eset <- imvigor210_eset[mad > 0.5, ]
#'pd1 <- as.numeric(imvigor210_eset["PDCD1", ])
#'group <- ifelse(pd1 > mean(pd1), "high", "low")
#'pd1_cor <- feature_select(x = imvigor210_eset, y = pd1, method = "cor", family = "pearson", padjcut = 0.05, cutoff = 0.5)
#'pd1_dif <- feature_select(x = imvigor210_eset, y = pd1, method = "dif", padjcut = 0.05, cutoff = 2)
#'pd1_dif_2 <- feature_select(x = imvigor210_eset, y = group, method = "dif", padjcut = 0.05, cutoff = 2)

feature_select <- function(x, y, method = c("cor", "dif"),
                           family = c("spearman", "pearson"),
                           cutoff = NULL, padjcut = NULL){
  method = match.arg(method)
  if (length(unique(y)) == 1){
    stop("There are only one constant value in y, y must be binary or quantitative value.")
  }
  type = ifelse(length(unique(y)) == 2, "Binary", "quantitative")
  if (length(unique(y)) == 2){
    message("Deteching two levels in y, we will treat y as a binary varibale")
  }
  if (length(unique(y)) > 2){
    message("Deteching more than two levels in y, we will treat y as a quantitative varibale")
  }

  if (method == "cor"){
    if (type != "quantitative"){
      stop("Correlation between x and y, y must be quantitative")
    }
    Gene = rownames(x)
    x <- as.tibble(t(x))
    tmp <- x %>%
      map(cor.test, y, method = family)
    pvalue <- tmp %>% map_dbl("p.value")
    estimate <- tmp %>% map_dbl("estimate")
    P.adj <- p.adjust(as.numeric(pvalue), method = "fdr")
    cor_feature <- tibble(Gene = Gene, P.adj = P.adj, Estimate = estimate)
    cor_feature <- cor_feature %>% filter(P.adj < padjcut, abs(Estimate) > cutoff) %>% select("Gene")
    feature <- cor_feature$Gene %>% as.character()
  }
  if (method == "dif"){
    if (type == "quantitative"){
      message("For quantitative varibale, upper 25% and bottom 25% samples
              were treated as upregulated group and downregulated group.")
      up <- which(y > quantile(y, 0.75)); down <- which(y < quantile(y, 0.25))
      pdata <- data.frame(samples = c(colnames(x)[up], colnames(x)[down]),
                          group = c(rep("up", length(up)), rep("down", length(down))))
      exprdata <- x[, c(up, down)]
      contrastfml <- c("up - down")
    }
    if (type == "Binary"){
      if (length(unique(y)) != 2){
        stop("y must be binary feature, please check your data carefully")
      }
      pdata <- data.frame(samples = colnames(x), group = y)
      contrastfml <- paste(unique(y)[1], "-", unique(y)[2])
      exprdata <- x
    }
    dif <- limma.dif(exprdata = exprdata, pdata = pdata, contrastfml = contrastfml)
    dif <- data.frame(Probe = rownames(dif), dif)
    dif_feature <- dif %>% as_tibble() %>% filter(abs(logFC) > cutoff, adj.P.Val < padjcut) %>%
      select("Probe")
    feature <- dif_feature$Probe %>% as.character()

  }
  return(feature)
}

#' limma.dif
#'
#' This function performs differential expression analysis using the limma package on a given gene expression dataset. 
#' It constructs a design matrix from phenotype data, fits a linear model, applies contrasts, and computes statistics for 
#' differential expression.
#'
#' @param exprdata input matrix.Rownames should be features like gene symbols or cgi, colnames be samples
#' @param pdata phenotype data.Two-column dataframe which column 1 should be the same with the colnames of exprdata and column 2 are the grouping variable.
#' @param contrastfml see ?makeContrasts
#' @import limma
#'
#' @return a data frame return by limma::toptable, which genes at rows and following columns: genelist, logFC, AveExpr...
#' @export
#'
#' @examples
limma.dif <- function(exprdata, pdata, contrastfml){
  group_list <- as.character(pdata[, 2])
  design <- model.matrix(~0 + factor(group_list))
  colnames(design) <- levels(as.factor(pdata[, 2]))
  rownames(design) <- colnames(exprdata)
  if (!all(colnames(exprdata) == pdata[, 1])){
    stop(" expression data do not match pdata")
  }
  contrast.matrix <- makeContrasts(contrasts = contrastfml, levels=design)
  fit <- lmFit(exprdata, design)
  fit <- contrasts.fit(fit,contrast.matrix)
  fit <- eBayes(fit)
  dif <- topTable(fit, adjust.method="BH", coef = contrastfml, number = Inf)
  return(dif)
}
