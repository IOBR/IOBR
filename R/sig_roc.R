#' Plot ROC Curves and Compare Them
#'
#' This function generates Receiver Operating Characteristic (ROC) curves for multiple predictors
#' and optionally performs statistical comparisons between them.
#'
#' @param data A data frame containing the predictor variables and the binary outcome variable.
#' @param response The name of the binary outcome variable in `data`.
#' @param variables A vector of names of predictor variables in `data` for which ROC curves will be plotted.
#' @param fig.path Directory path to save the output PDF file. Default is the current working directory.
#' @param main Main title for the ROC plot.
#' @param file.name Name of the output PDF file without extension. Default is "0-ROC of multiple variables".
#' @param palette Color palette for plotting ROC curves. Default is "jama".
#' @param cols Optional vector of colors for ROC curves. If NULL, colors are assigned automatically.
#' @param alpha Transparency level of colors (1 = opaque, 0 = transparent). Default is 1.
#' @param compare Logical indicating whether to perform statistical comparison of AUCs. Default is FALSE.
#' @param smooth Logical indicating whether to smooth ROC curves. Default is TRUE.
#' @param compare_method Method for comparing ROC curves if `compare` is TRUE. Default is "bootstrap".
#' @param boot.n Number of bootstrap replications for comparison. Default is 100.
#'
#' @return A list containing:
#'   - `auc.out`: Data frame with AUC values and confidence intervals for each variable.
#'   - `legend.name`: Vector of legend entries for the plot.
#'   - `p.out`: If `compare` is TRUE, data frame with p-values from pairwise comparisons.
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' data("tcga_stad_pdata", package = "IOBR")
#' sig_roc(data = tcga_stad_pdata, response = "OS_status",
#'         variables = c("TMEscore_plus", "GZMB", "GNLY"))
sig_roc <- function(data,
                    response,
                    variables,
                    fig.path = ".",
                    main = NULL,
                    file.name = NULL,
                    palette = "jama",
                    cols = NULL,
                    alpha = 1,
                    compare = FALSE,
                    smooth = TRUE,
                    compare_method = "bootstrap",
                    boot.n = 100) {
  # Store current pROC options and set progress to none
  old_options <- getOption("pROCProgress")
  on.exit(options(pROCProgress = old_options))
  options(pROCProgress = list(name = "none"))


  if (!response %in% colnames(data)) {
    stop("The response column specified does not exist in the data.")
  }

  data <- as.data.frame(data)

  # Ensure response is a factor with two levels
  data[[response]] <- as.factor(data[[response]])
  if (length(levels(data[[response]])) != 2) {
    stop("Response variable must have exactly two levels.")
  }

  data <- data[!is.na(data[[response]]), ]
  variables <- variables[variables %in% colnames(data)]


  input <- as.data.frame(data[, c(response, variables)])

  message(">>>== head of input data: ")
  print(head(input))

  var_counts <- length(variables[variables %in% colnames(data)])

  if (is.null(cols)) {
    cols <- palettes(palette = palette, alpha = alpha, show_col = FALSE, show_message = FALSE)

    if (var_counts > length(cols)) {
      cols <- palettes(category = "random", alpha = alpha, show_col = FALSE, show_message = FALSE)
    }
  }

  ########################################

  auc.out <- c()

  if (is.null(file.name)) file.name <- "0-ROC of multiple variables"

  outfile <- file.path(fig.path, paste(file.name, ".pdf", sep = ""))

  pdf(file = outfile, width = 5, height = 5)

  on.exit(dev.off(), add = TRUE) # Ensure graphics device is shut off when function exits
  x <- pROC::plot.roc(input[, 1], input[, 2],
    ylim = c(0, 1),
    xlim = c(1, 0),
    smooth = smooth,
    ci = TRUE,
    main = main,
    col = cols[2],
    lwd = 1.5,
    legacy.axes = T,
    xlab = "False Positive Rate",
    ylab = "True Positive Rate"
  )

  ci.lower <- round(as.numeric(x$ci[1]), 3)
  ci.upper <- round(as.numeric(x$ci[3]), 3)

  auc.ci <- c(colnames(input)[2], round(as.numeric(x$auc), 3), paste(ci.lower, ci.upper, sep = "-"))
  auc.out <- rbind(auc.out, auc.ci)
  ##############################
  for (i in 3:ncol(input)) {
    x <- pROC::plot.roc(input[, 1], input[, i],
      add = T,
      smooth = smooth,
      ci = TRUE,
      col = cols[i],
      lwd = 2,
      legacy.axes = T,
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

  legend.name <- paste(colnames(input)[2:length(input)], " AUC = ", auc.out$AUC, sep = " ")
  legend("bottomright",
    legend = legend.name,
    col = cols[2:length(input)],
    lwd = 2,
    bty = "n"
  )

  # dev.off()
  ###################################################
  if (compare) {
    p.out <- c()

    for (i in 2:(ncol(input) - 1)) {
      for (j in (i + 1):ncol(input)) {
        p <- pROC::roc.test(input[, 1], input[, i], input[, j], method = compare_method, boot.n = boot.n, progress = "none")
        p.tmp <- c(colnames(input)[i], colnames(input)[j], p$p.value)
        p.out <- rbind(p.out, p.tmp)
      }
    }
    p.out <- as.data.frame(p.out)
    # head(p.out)
    colnames(p.out) <- c("ROC1", "ROC2", "p.value")
    p.out$p.value <- round(as.numeric(p.out$p.value), 5)

    return(list(auc.out = auc.out, legend.name = legend.name, p.out = p.out))
  } else {
    return(list(auc.out = auc.out, legend.name = legend.name))
  }
}
