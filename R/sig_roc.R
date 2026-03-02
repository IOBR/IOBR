#' Plot ROC Curves and Compare Them
#'
#' This function generates Receiver Operating Characteristic (ROC) curves for multiple predictors
#' and optionally performs statistical comparisons between them.
#'
#' @param data A data frame containing the predictor variables and the binary outcome variable.
#' @param response The name of the binary outcome variable in `data`.
#' @param variables A vector of names of predictor variables in `data` for which ROC curves will be plotted.
#' @param fig.path Directory path to save the output PDF file. Default is NULL.
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

  # var_counts <- length(variables[variables %in% colnames(data)])
  
  # 数据准备
  var_names <- variables  # 明确变量名（已过滤存在的列）
  n_vars <- length(var_names)

  if (is.null(cols)) {
    cols <- palettes(palette = palette, alpha = alpha, show_col = FALSE, show_message = FALSE)

    if (n_vars > length(cols)) {
      cols <- palettes(category = "random", alpha = alpha, show_col = FALSE, show_message = FALSE)
    }
  }

  ########################################

  auc.out <- c()

  # if (is.null(file.name)) file.name <- "0-ROC of multiple variables"
  # 
  # outfile <- file.path(fig.path, paste(file.name, ".pdf", sep = ""))
  # 
  # pdf(file = outfile, width = 5, height = 5)
  # 
  # on.exit(dev.off(), add = TRUE) # Ensure graphics device is shut off when function exits

  rlang::check_installed("pROC")
  # x <- pROC::plot.roc(input[, 1], input[, 2],
  #   ylim = c(0, 1),
  #   xlim = c(1, 0),
  #   smooth = smooth,
  #   ci = TRUE,
  #   main = main,
  #   col = cols[2],
  #   lwd = 1.5,
  #   legacy.axes = T,
  #   xlab = "False Positive Rate",
  #   ylab = "True Positive Rate"
  # )
  # 
  # ci.lower <- round(as.numeric(x$ci[1]), 3)
  # ci.upper <- round(as.numeric(x$ci[3]), 3)
  # 
  # auc.ci <- c(colnames(input)[2], round(as.numeric(x$auc), 3), paste(ci.lower, ci.upper, sep = "-"))
  # auc.out <- rbind(auc.out, auc.ci)
  # ##############################
  # for (i in 3:ncol(input)) {
  #   x <- pROC::plot.roc(input[, 1], input[, i],
  #     add = T,
  #     smooth = smooth,
  #     ci = TRUE,
  #     col = cols[i],
  #     lwd = 2,
  #     legacy.axes = T,
  #     xlab = "False Positive Rate",
  #     ylab = "True Positive Rate"
  #   )
  # 
  #   ci.lower <- round(as.numeric(x$ci[1]), 3)
  #   ci.upper <- round(as.numeric(x$ci[3]), 3)
  #   auc.ci <- c(colnames(input)[i], round(as.numeric(x$auc), 3), paste(ci.lower, ci.upper, sep = "-"))
  #   auc.out <- rbind(auc.out, auc.ci)
  # }
  # 
  # auc.out <- as.data.frame(auc.out)
  # colnames(auc.out) <- c("Name", "AUC", "AUC CI")
  # 
  # legend.name <- paste(colnames(input)[2:length(input)], " AUC = ", auc.out$AUC, sep = " ")
  # legend("bottomright",
  #   legend = legend.name,
  #   col = cols[2:length(input)],
  #   lwd = 2,
  #   bty = "n"
  # )

  # dev.off()
  

  cols_use <- cols[1:n_vars]  # 明确取前 n 个
  
  # 新增：收集 ROC 数据用于 ggplot2
  roc_objects <- list()
  roc_data_list <- list()
  auc_data <- data.frame(
    Name = character(),
    AUC = numeric(),
    AUC_CI = character(),
    stringsAsFactors = FALSE
  )
  
  # 计算所有 ROC 对象
  for (i in seq_along(var_names)) {
    var_name <- var_names[i]
    
    roc_obj <- pROC::roc(input[[response]], input[[var_name]], 
                         smooth = smooth, ci = TRUE)
    roc_objects[[var_name]] <- roc_obj
    
    # 提取曲线数据（注意：specificity 是 1-FPR，需要转换）
    roc_df <- data.frame(
      FPR = 1 - roc_obj$specificities,   # False Positive Rate
      TPR = roc_obj$sensitivities,       # True Positive Rate
      variable = var_name,
      stringsAsFactors = FALSE
    )
    roc_data_list[[i]] <- roc_df
    
    # 提取 AUC
    ci_vals <- as.numeric(roc_obj$ci)
    auc_data <- rbind(auc_data, data.frame(
      Name = var_name,
      AUC = round(as.numeric(roc_obj$auc), 3),
      AUC_CI = paste(round(ci_vals[1], 3), round(ci_vals[3], 3), sep = "-")
    ))
  }
  
  # 合并数据
  all_roc_data <- do.call(rbind, roc_data_list)
  
  # 确保变量顺序与颜色对应
  all_roc_data$variable <- factor(all_roc_data$variable, levels = var_names)
  
  ########################################
  
  
  p <- ggplot2::ggplot(all_roc_data, 
                       ggplot2::aes(x = FPR, y = TPR, color = variable)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, 
                         linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(
      values = cols_use,
      labels = paste0(auc_data$Name, " (AUC=", auc_data$AUC, " [", auc_data$AUC_CI, "])"),
      name = "Variables"
    ) +
    ggplot2::labs(
      x = "False Positive Rate",
      y = "True Positive Rate",
      title = main
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottomright",
      legend.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  ########################################
  # 保存
  if (!is.null(fig.path)) {
    if (is.null(file.name)) {
      file.name <- "0-ROC of multiple variables"
    }
    if (!dir.exists(fig.path)) {
      dir.create(fig.path, recursive = TRUE, showWarnings = FALSE)
    }
    outfile <- file.path(fig.path, paste0(file.name, ".pdf"))
    ggplot2::ggsave(outfile, plot = p, width = 5, height = 5)
    message(">>>-- Plot saved to: ", outfile)
  }
  
  ########################################
  # 返回
  legend.name <- paste0(auc_data$Name, " AUC = ", auc_data$AUC)
  auc.out <- auc_data
  colnames(auc.out) <- c("Name", "AUC", "AUC CI")
  
  result <- list(
    plot = p,
    auc.out = auc.out,
    legend.name = legend.name
  )
  
  # 比较
  if (compare) {
    p.out <- data.frame(
      ROC1 = character(),
      ROC2 = character(),
      p.value = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:(n_vars - 1)) {
      for (j in (i + 1):n_vars) {
        test_res <- pROC::roc.test(
          roc_objects[[var_names[i]]], 
          roc_objects[[var_names[j]]],
          method = compare_method,
          boot.n = boot.n,
          progress = "none"
        )
        p.out <- rbind(p.out, data.frame(
          ROC1 = var_names[i],
          ROC2 = var_names[j],
          p.value = round(test_res$p.value, 5),
          stringsAsFactors = FALSE
        ))
      }
    }
    result$p.out <- p.out
  }
  
  return(result)
  
  ###################################################
  # if (compare) {
  #   p.out <- c()
  # 
  #   for (i in 2:(ncol(input) - 1)) {
  #     for (j in (i + 1):ncol(input)) {
  #       p <- pROC::roc.test(input[, 1], input[, i], input[, j], method = compare_method, boot.n = boot.n, progress = "none")
  #       p.tmp <- c(colnames(input)[i], colnames(input)[j], p$p.value)
  #       p.out <- rbind(p.out, p.tmp)
  #     }
  #   }
  #   p.out <- as.data.frame(p.out)
  #   # head(p.out)
  #   colnames(p.out) <- c("ROC1", "ROC2", "p.value")
  #   p.out$p.value <- round(as.numeric(p.out$p.value), 5)
  # 
  #   return(list(auc.out = auc.out, legend.name = legend.name, p.out = p.out))
  # } else {
  #   return(list(auc.out = auc.out, legend.name = legend.name))
  # }
  
}
