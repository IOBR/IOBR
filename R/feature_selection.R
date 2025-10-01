#' Feature Selection via Correlation or Differential Expression
#'
#' Selects informative features using either correlation with a quantitative response or differential expression (limma) for binary/continuous responses.
#'
#' @param x Numeric matrix. Features (rows) by samples (columns).
#' @param y Numeric or factor. Response vector (quantitative or binary).
#' @param method Character. "cor" (correlation) or "dif" (differential expression). Default c("cor","dif").
#' @param family Character. Correlation method if method = "cor": "spearman" or "pearson".
#' @param cutoff Numeric. Absolute correlation (for cor) or |log2FC| (for dif) threshold.
#' @param padjcut Numeric. Adjusted p-value cutoff.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import purrr
#' @import stringr
#' @import glmnet
#'
#' @return Character vector of selected feature names.
#' @export
#'
#' @examples
#' data("imvigor210_eset", package = "IOBR")
#' mad <- apply(imvigor210_eset, 1, mad)
#' imvigor210_eset <- imvigor210_eset[mad > 0.5, ]
#' pd1 <- as.numeric(imvigor210_eset["PDCD1", ])
#' group <- ifelse(pd1 > mean(pd1), "high", "low")
#' pd1_cor <- feature_select(x = imvigor210_eset, y = pd1, method = "cor", family = "pearson", padjcut = 0.05, cutoff = 0.5)
#' pd1_dif <- feature_select(x = imvigor210_eset, y = pd1, method = "dif", padjcut = 0.05, cutoff = 2)
#' pd1_dif_2 <- feature_select(x = imvigor210_eset, y = group, method = "dif", padjcut = 0.05, cutoff = 2)
feature_select <- function(x, y, method = c("cor", "dif"),
                           family = c("spearman", "pearson"),
                           cutoff = NULL, padjcut = NULL) {
  method <- match.arg(method)
  if (length(unique(y)) == 1) {
    stop("There are only one constant value in y, y must be binary or quantitative value.")
  }
  type <- ifelse(length(unique(y)) == 2, "Binary", "quantitative")
  if (length(unique(y)) == 2) {
    message("Deteching two levels in y, we will treat y as a binary varibale")
  }
  if (length(unique(y)) > 2) {
    message("Deteching more than two levels in y, we will treat y as a quantitative varibale")
  }

  if (method == "cor") {
    if (type != "quantitative") {
      stop("Correlation between x and y, y must be quantitative")
    }
    Gene <- rownames(x)
    x <- as.tibble(t(x))
    tmp <- x %>%
      map(cor.test, y, method = family)
    pvalue <- tmp %>% map_dbl("p.value")
    estimate <- tmp %>% map_dbl("estimate")
    P.adj <- p.adjust(as.numeric(pvalue), method = "fdr")
    cor_feature <- tibble(Gene = Gene, P.adj = P.adj, Estimate = estimate)
    cor_feature <- cor_feature %>%
      filter(P.adj < padjcut, abs(Estimate) > cutoff) %>%
      select("Gene")
    feature <- cor_feature$Gene %>% as.character()
  }
  if (method == "dif") {
    if (type == "quantitative") {
      message("For quantitative varibale, upper 25% and bottom 25% samples
              were treated as upregulated group and downregulated group.")
      up <- which(y > quantile(y, 0.75))
      down <- which(y < quantile(y, 0.25))
      pdata <- data.frame(
        samples = c(colnames(x)[up], colnames(x)[down]),
        group = c(rep("up", length(up)), rep("down", length(down)))
      )
      exprdata <- x[, c(up, down)]
      contrastfml <- c("up - down")
    }
    if (type == "Binary") {
      if (length(unique(y)) != 2) {
        stop("y must be binary feature, please check your data carefully")
      }
      pdata <- data.frame(samples = colnames(x), group = y)
      contrastfml <- paste(unique(y)[1], "-", unique(y)[2])
      exprdata <- x
    }
    dif <- limma.dif(exprdata = exprdata, pdata = pdata, contrastfml = contrastfml)
    dif <- data.frame(Probe = rownames(dif), dif)
    dif_feature <- dif %>%
      as_tibble() %>%
      filter(abs(logFC) > cutoff, adj.P.Val < padjcut) %>%
      select("Probe")
    feature <- dif_feature$Probe %>% as.character()
  }
  return(feature)
}

#' Differential Expression Analysis Using Limma
#'
#' Performs differential expression analysis using the limma package on a given gene expression dataset.
#' Constructs a design matrix from phenotype data, fits a linear model, applies contrasts, and computes
#' statistics for differential expression.
#'
#' @param exprdata A matrix with rownames as features like gene symbols or cgi, and colnames as samples.
#' @param pdata A two-column dataframe where the first column matches the colnames of exprdata and the second column contains the grouping variable.
#' @param contrastfml A character vector for contrasts to be tested (see ?makeContrasts for more details).
#'
#' @return Returns a dataframe from limma::topTable, which includes genes as rows and columns like genelist, logFC, AveExpr, etc.
#' @export
#' @examples
#' data("expression_data", package = "ExamplePackage")
#' data("phenotype_data", package = "ExamplePackage")
#'
#' dif_results <- limma.dif(exprdata = expression_data, pdata = phenotype_data, contrastfml = "group1 - group2")
#' print(dif_results)
limma.dif <- function(exprdata, pdata, contrastfml) {
  group_list <- as.character(pdata[, 2])
  design <- model.matrix(~ 0 + factor(group_list))
  colnames(design) <- levels(as.factor(pdata[, 2]))
  rownames(design) <- colnames(exprdata)
  if (!all(colnames(exprdata) == pdata[, 1])) {
    stop(" expression data do not match pdata")
  }
  contrast.matrix <- makeContrasts(contrasts = contrastfml, levels = design)
  fit <- lmFit(exprdata, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  dif <- topTable(fit, adjust.method = "BH", coef = contrastfml, number = Inf)
  return(dif)
}
