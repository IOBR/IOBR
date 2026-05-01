#' Plot ROC Curves and Compare Them
#'
#' @description
#' Generates Receiver Operating Characteristic (ROC) curves for multiple
#' predictors and optionally performs statistical comparisons between them.
#'
#' @param data Data frame containing the predictor variables and binary outcome.
#' @param response Character. Name of the binary outcome variable in `data`.
#' @param variables Character vector. Names of predictor variables for ROC curves.
#' @param fig.path Character or `NULL`. Directory path to save output PDF.
#'   Default is `NULL`.
#' @param main Character or `NULL`. Main title for the ROC plot. Default is `NULL`.
#' @param file.name Character or `NULL`. Output PDF filename without extension.
#'   Default is `"0-ROC of multiple variables"`.
#' @param palette Character. Color palette for ROC curves. Default is `"jama"`.
#' @param cols Character vector or `NULL`. Custom colors for ROC curves.
#'   Default is `NULL`.
#' @param alpha Numeric. Transparency level (1 = opaque, 0 = transparent).
#'   Default is `1`.
#' @param compare Logical. Whether to perform statistical comparison of AUCs.
#'   Default is `FALSE`.
#' @param smooth Logical. Whether to smooth ROC curves. Default is `TRUE`.
#' @param compare_method Character. Method for comparing ROC curves.
#'   Default is `"bootstrap"`.
#' @param boot.n Integer. Number of bootstrap replications. Default is `100`.
#'
#' @return A list containing:
#' \describe{
#'   \item{auc.out}{Data frame with AUC values and confidence intervals}
#'   \item{legend.name}{Vector of legend entries for the plot}
#'   \item{p.out}{If `compare = TRUE`, data frame with p-values from comparisons}
#' }
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' tcga_stad_pdata <- load_data("tcga_stad_pdata")
#' sig_roc(
#'   data = tcga_stad_pdata, response = "OS_status",
#'   variables = c("TMEscore_plus", "GZMB", "GNLY")
#' )
#' }
sig_roc <- function(data,
                    response,
                    variables,
                    fig.path = NULL,
                    main = NULL,
                    file.name = NULL,
                    palette = "jama",
                    cols = NULL,
                    alpha = 1,
                    compare = FALSE,
                    smooth = TRUE,
                    compare_method = "bootstrap",
                    boot.n = 100) {
  options(pROCProgress = list(name = "none"))

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  if (!response %in% colnames(data)) {
    cli::cli_abort("Response column {.val {response}} not found in data.")
  }

  data <- as.data.frame(data)
  data[[response]] <- as.factor(data[[response]])

  if (length(levels(data[[response]])) != 2) {
    cli::cli_abort("Response variable must have exactly two levels.")
  }

  data <- data[!is.na(data[[response]]), ]
  variables <- variables[variables %in% colnames(data)]

  if (length(variables) == 0) {
    cli::cli_abort("No valid variables found in data.")
  }

  input <- as.data.frame(data[, c(response, variables)])

  cli::cli_alert_info("Input data preview:")
  if (interactive()) print(head(input))

  var_counts <- length(variables)

  if (is.null(cols)) {
    cols <- palettes(palette = palette, alpha = alpha, show_col = FALSE, show_message = FALSE)
    if (var_counts > length(cols)) {
      cols <- palettes(category = "random", alpha = alpha, show_col = FALSE, show_message = FALSE)
    }
  }

  auc.out <- c()

  if (!is.null(fig.path)) {
    if (is.null(file.name)) file.name <- "0-ROC of multiple variables"
    outfile <- file.path(fig.path, paste0(file.name, ".pdf"))
    if (!dir.exists(fig.path)) {
      dir.create(fig.path, recursive = TRUE, showWarnings = FALSE)
    }
    pdf(file = outfile, width = 5, height = 5)
    on.exit(dev.off(), add = FALSE)
  }

  x <- pROC::plot.roc(input[, 1], input[, 2],
    #ylim = c(0, 1),
    #xlim = c(1, 0),
    smooth = smooth,
    ci = TRUE,
    main = main,
    col = cols[2],
    lwd = 1.5,
    legacy.axes = TRUE,
    xlab = "False Positive Rate",
    ylab = "True Positive Rate"
  )

  ci.lower <- round(as.numeric(x$ci[1]), 3)
  ci.upper <- round(as.numeric(x$ci[3]), 3)

  auc.ci <- c(colnames(input)[2], round(as.numeric(x$auc), 3), paste(ci.lower, ci.upper, sep = "-"))
  auc.out <- rbind(auc.out, auc.ci)

  for (i in seq(3, ncol(input))) {
    x <- pROC::plot.roc(input[, 1], input[, i],
      add = TRUE,
      smooth = smooth,
      ci = TRUE,
      col = cols[i],
      lwd = 2,
      legacy.axes = TRUE,
      xlab = "False Positive Rate",
      ylab = "True Positive Rate"
    )

    ci.lower <- round(as.numeric(x$ci[1]), 3)
    ci.upper <- round(as.numeric(x$ci[3]), 3)
    auc.ci <- c(colnames(input)[i], round(as.numeric(x$auc), 3), paste(ci.lower, ci.upper, sep = "-"))
    auc.out <- rbind(auc.out, auc.ci)
  }

  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name", "AUC", "AUC CI")

  legend.name <- paste(colnames(input)[seq(2, ncol(input))], " AUC = ", auc.out$AUC, sep = " ")
  #oldpar <- par(xpd = TRUE)
  #on.exit(par(oldpar), add = TRUE)
  legend("bottomright",
    legend = legend.name,
    col = cols[seq(2, length(variables) + 1)],
    lwd = 2,
    bty = "n"
  )

  if (compare) {
    p.out <- c()
    for (i in seq(2, ncol(input) - 1)) {
      for (j in seq(i + 1, ncol(input))) {
        p <- pROC::roc.test(input[, 1], input[, i], input[, j],
          method = compare_method, boot.n = boot.n, progress = "none"
        )
        p.tmp <- c(colnames(input)[i], colnames(input)[j], p$p.value)
        p.out <- rbind(p.out, p.tmp)
      }
    }
    p.out <- as.data.frame(p.out)
    colnames(p.out) <- c("ROC1", "ROC2", "p.value")
    p.out$p.value <- round(as.numeric(p.out$p.value), 5)
    return(list(auc.out = auc.out, legend.name = legend.name, p.out = p.out))
  } else {
    return(list(auc.out = auc.out, legend.name = legend.name))
  }
}
