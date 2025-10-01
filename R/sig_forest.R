#' Forest Plot for Survival Analysis Results
#'
#' Generates a forest plot to visualize hazard ratios, confidence intervals, and p-values
#' for gene signatures or features from survival analysis.
#'
#' @param data Data frame with survival analysis results including p-values, HRs, and CIs.
#' @param signature Column name for signatures or feature names.
#' @param pvalue Column name for p-values. Default is "P".
#' @param HR Column name for hazard ratios. Default is "HR".
#' @param CI_low_0.95 Column name for lower CI bound. Default is "CI_low_0.95".
#' @param CI_up_0.95 Column name for upper CI bound. Default is "CI_up_0.95".
#' @param n Maximum number of signatures to display. Default is 10.
#' @param max_character Maximum characters for labels before wrapping. Default is 25.
#' @param discrete_width Width for discretizing long labels. Default is 35.
#' @param color_option Color option for p-value gradient (1, 2, or 3). Default is 1.
#' @param text.size Text size for y-axis labels. Default is 13.
#'
#' @return A ggplot2 object of the forest plot.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' sig_surv_result <- batch_surv(pdata = pdata_sig_tme_binary, variable = c(100:ncol(pdata_sig_tme_binary)))
#' sig_forest(data = sig_surv_result, signature = "ID")
sig_forest <- function(data, signature, pvalue = "P", HR = "HR", CI_low_0.95 = "CI_low_0.95",
                       CI_up_0.95 = "CI_up_0.95", n = 10, max_character = 25,
                       discrete_width = 35, color_option = 1,
                       text.size = 13) {
  data <- as.data.frame(data)
  colnames(data)[which(colnames(data) == signature)] <- "signature"
  colnames(data)[which(colnames(data) == pvalue)] <- "P"
  colnames(data)[which(colnames(data) == HR)] <- "HR"
  colnames(data)[which(colnames(data) == CI_low_0.95)] <- "CI_low_0.95"
  colnames(data)[which(colnames(data) == CI_up_0.95)] <- "CI_up_0.95"

  data[, c("P", "HR", "CI_low_0.95", "CI_up_0.95")] <- apply(data[, c("P", "HR", "CI_low_0.95", "CI_up_0.95")], 2, as.numeric)
  data <- data[complete.cases(data), ]

  if (dim(data)[1] > n) {
    message(paste0("Top ", n, " signatures will be shown"))

    data$HR_statistic <- data$HR - 1
    good_features <- high_var_fea(
      result = data,
      target = "signature",
      name_padj = "P",
      padj_cutoff = 1,
      name_logfc = "HR_statistic",
      logfc_cutoff = 0,
      n = n / 2
    )
    data <- data[data$signature %in% good_features, ]
    data <- data[order(data$HR, decreasing = FALSE), ]
    data$signature <- as.character(data$signature)
  }

  if (max(nchar(data$signature)) > max_character) {
    goi <- as.character(data$signature)
    for (i in 1:length(goi)) {
      if (nchar(goi[i]) > max_character) {
        data[i, "signature"] <- gsub(data[i, "signature"], pattern = "\\_", replacement = " ")
      }
    }
  }

  # Set the order of 'signature' as a factor based on 'HR'
  data$signature <- factor(data$signature, levels = data$signature)

  pp <- ggplot(data = data, aes(x = HR, y = signature, color = P)) +
    geom_errorbarh(aes(xmax = CI_up_0.95, xmin = CI_low_0.95), color = "black", height = 0, size = 1.2) +
    geom_point(aes(x = HR, y = signature), size = 4.5, shape = 16) +
    geom_vline(xintercept = 1, linetype = "dashed", size = 0.8) +
    scale_x_continuous(breaks = c(0.5, 1, 1.50)) +
    coord_trans(x = "log2") +
    ylab("Features") +
    xlab("Hazard Ratios of Features") +
    labs(color = "P value") +
    viridis::scale_color_viridis(option = color_option) +
    theme_light() +
    theme(
      axis.text.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0),
      axis.text.y = element_text(size = text.size, color = "black", vjust = 0.5, hjust = 1, angle = 0),
      title = element_text(size = 15, colour = "black", vjust = 0.5, hjust = 0.5)
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = discrete_width))

  return(pp)
}
