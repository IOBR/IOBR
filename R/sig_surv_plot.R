#' Generate Kaplan-Meier Survival Plot for Signature
#'
#' @description
#' Creates Kaplan-Meier survival plots for a given signature or gene, with
#' automatic cutoff determination. Generates three types of plots: optimal
#' cutoff (best cutoff), tertile-based (3 groups), and median split (2 groups).
#'
#' @param input_pdata Data frame with survival data and signature scores.
#' @param signature Character string. Column name of the target signature.
#' @param project Character string. Project name for output. Default is `"KM"`.
#' @param time Character string. Column name for survival time. Default is `"time"`.
#' @param status Character string. Column name for survival status. Default is `"status"`.
#' @param time_type Character string. Time unit (`"month"` or `"day"`).
#'   Default is `"month"`.
#' @param break_month Numeric or `"auto"`. Time axis breaks. Default is `"auto"`.
#' @param palette Character string. Color palette if `cols` not provided.
#'   Default is `"jama"`.
#' @param save_path Character string or `NULL`. Directory for saving plots.
#'   If `NULL`, plots are not saved. Default is `NULL`.
#' @param mini_sig Character string. Label for low score group. Default is `"score"`.
#' @param show_col Logical indicating whether to show colors. Default is `TRUE`.
#' @param index Integer. Index for multiple plots. Default is `1`.
#' @param cols Character vector. Optional custom colors.
#' @param fig.type Character string. File format. Default is `"png"`.
#' @param ID Character string. Column name for sample IDs. Default is `"ID"`.
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{Processed input data with group assignments}
#'   \item{plots}{Combined survival plots}
#' }
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' tcga_stad_pdata <- load_data("tcga_stad_pdata")
#' sig_surv_plot(
#'   input_pdata = tcga_stad_pdata,
#'   signature = "TMEscore_plus",
#'   time = "time",
#'   status = "OS_status"
#' )
sig_surv_plot <- function(input_pdata,
                          signature,
                          project = "KM",
                          ID = "ID",
                          time = "time",
                          status = "status",
                          time_type = "month",
                          break_month = "auto",
                          cols = NULL,
                          palette = "jama",
                          show_col = TRUE,
                          mini_sig = "score",
                          fig.type = "png",
                          save_path = NULL,
                          index = 1) {
  rlang::check_installed("survminer")

  # Input validation
  if (!signature %in% colnames(input_pdata)) {
    cli::cli_abort("Signature {.val {signature}} not found in input_pdata")
  }

  # Create output directory
  abspath <- NULL
  if (!is.null(save_path)) {
    cli::cli_alert_info("Saving plots to: {.val {save_path}}")
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    abspath <- file.path(normalizePath(save_path, winslash = "/", mustWork = FALSE), "")
  }

  # Select relevant columns
  input_pdata <- as.data.frame(input_pdata[, colnames(input_pdata) %in% c(ID, time, status, signature)])

  # Standardize column names
  colnames(input_pdata)[colnames(input_pdata) == time] <- "time"
  colnames(input_pdata)[colnames(input_pdata) == status] <- "status"
  colnames(input_pdata)[colnames(input_pdata) == ID] <- "ID"

  # Convert and filter data
  input_pdata$time <- as.numeric(input_pdata$time)
  input_pdata$status <- as.numeric(as.character(input_pdata$status))

  input_pdata <- input_pdata %>%
    dplyr::filter(!is.na(.data$time)) %>%
    dplyr::filter(!is.na(.data$status)) %>%
    dplyr::filter(.data$time > 0)

  input_pdata$ID <- as.character(input_pdata$ID)

  # Transform follow-up time to months if needed
  if (time_type == "day") {
    input_pdata$time <- input_pdata$time / 30
  }

  # Create tertile-based groups
  q1 <- stats::quantile(input_pdata[[signature]], probs = 1 / 3, na.rm = TRUE)
  q2 <- stats::quantile(input_pdata[[signature]], probs = 2 / 3, na.rm = TRUE)

  input_pdata$group3 <- ifelse(
    input_pdata[[signature]] <= q1, "Low",
    ifelse(input_pdata[[signature]] >= q2, "High", "Middle")
  )

  # Create mean-based groups
  mean_val <- mean(input_pdata[[signature]], na.rm = TRUE)
  input_pdata$group2 <- ifelse(input_pdata[[signature]] <= mean_val, "Low", "High")

  cli::cli_alert_info(
    "Survival follow-up time range: {round(min(input_pdata$time, na.rm = TRUE), 2)} to {round(max(input_pdata$time, na.rm = TRUE), 2)} months"
  )

  # Calculate best cutoff
  res.cut <- survminer::surv_cutpoint(
    input_pdata,
    time = "time", event = "status", variables = signature
  )
  res.cut <- res.cut$cutpoint[[1]]
  cli::cli_alert_info("Best cutoff for {.val {signature}}: {round(res.cut, 2)}")

  # Apply best cutoff
  input_pdata <- best_cutoff(
    pdata = input_pdata, variable = signature,
    time = "time", status = "status", print_result = FALSE
  )

  colnames(input_pdata)[colnames(input_pdata) == paste0(signature, "_binary")] <- "bestcutoff"

  cli::cli_alert_info("High {signature}: {table(input_pdata$bestcutoff)['High']}")
  cli::cli_alert_info("Low {signature}: {table(input_pdata$bestcutoff)['Low']}")

  # Save data if path provided
  if (!is.null(save_path)) {
    save(input_pdata, file = file.path(
      save_path,
      paste0(index, "-0-", project, "-", signature, "-survival-analysis-input.RData")
    ))
  }

  # Prepare data for Cox model (High = 1)
  input_pdata$bestcutoff <- ifelse(input_pdata$bestcutoff == "High", 1, 0)

  # Calculate HR and CI
  cox_fit <- survival::coxph(
    survival::Surv(time = input_pdata$time, event = input_pdata$status) ~ bestcutoff,
    data = input_pdata
  )
  pvalue <- getHRandCIfromCoxph(cox_fit)

  HR <- paste0("Hazard Ratio = ", round(pvalue[, 2], 2))
  CI <- paste0("95% CI: ", round(pvalue[, 3], 2), " - ", round(pvalue[, 4], 2))
  cut_off <- paste0("cutoff = ", round(res.cut, 3))

  # Create survival fit
  sfit <- survminer::surv_fit(
    survival::Surv(time = input_pdata$time, event = input_pdata$status) ~ bestcutoff,
    data = input_pdata
  )

  # Define break time
  break_month_val <- break_month
  if (identical(break_month_val, "auto")) {
    break_month_val <- calculate_break_month(input = input_pdata$time, block = 6)
  }
  break_month_val <- as.numeric(break_month_val)
  max_month <- break_month_val * 6

  # Get colors
  if (is.null(cols)) {
    cols <- palettes(category = "box", palette = palette, show_col = FALSE)
  }

  # Plot 1: Best cutoff
  pp1 <- survminer::ggsurvplot(
    sfit,
    data = input_pdata, censor = TRUE, ncensor.plot = FALSE, conf.int = FALSE,
    xlim = c(0, max_month), break.time.by = break_month_val,
    xlab = "Months after diagnosis", surv.median.line = "h",
    legend.labs = c(paste0("Low ", mini_sig), paste0("High ", mini_sig)),
    submain = paste0(signature, "-in-", project),
    risk.table = TRUE, tables.height = 0.20, palette = cols, pval.size = 8,
    pval = paste(
      ifelse(pvalue[, 1] < 0.0001, "P < 0.0001",
        paste0("P = ", round(pvalue[, 1], 4))
      ),
      HR, CI, cut_off,
      sep = "\n"
    ),
    size = 0.4
  )

  res1 <- survminer::arrange_ggsurvplots(list(pp1), print = FALSE, ncol = 1, nrow = 1)

  if (!is.null(save_path)) {
    ggplot2::ggsave(
      plot = res1,
      filename = paste0(
        index, "-1-KMplot-best-cutoff-", signature, "-",
        project, ".", fig.type
      ),
      width = 6, height = 6.5, path = save_path
    )
  }

  # Plot 2: Three groups (tertiles)
  input_pdata$bestcutoff <- ifelse(input_pdata$bestcutoff == 1, "High", "Low")

  input_pdata <- input_pdata[!is.na(input_pdata$group3), ]
  sfit <- survminer::surv_fit(
    survival::Surv(input_pdata$time, input_pdata$status) ~ group3,
    data = input_pdata
  )

  names(sfit$strata) <- gsub("group3=", "", names(sfit$strata))

  pp2 <- survminer::ggsurvplot(
    sfit,
    data = input_pdata, censor = TRUE, ncensor.plot = FALSE,
    conf.int = FALSE, xlim = c(0, max_month), break.time.by = break_month_val,
    xlab = "Months after diagnosis",
    submain = paste0(signature, "-in-", project),
    surv.median.line = "h", risk.table = TRUE, tables.height = 0.25,
    palette = cols, pval.size = 8, size = 0.4
  )

  # Log-rank test
  fitd <- survival::survdiff(
    survival::Surv(time, status) ~ group3,
    data = input_pdata,
    na.action = stats::na.exclude
  )
  p.val <- 1 - stats::pchisq(fitd$chisq, length(fitd$n) - 1)

  p.lab <- paste0(
    "Overall P",
    ifelse(p.val < 0.001, " < 0.001", paste0(" = ", round(p.val, 3)))
  )

  pp2$plot <- pp2$plot + ggplot2::annotate(
    "text",
    x = 0, y = 0.55, hjust = 0, fontface = 3, label = p.lab
  )

  # Pairwise comparisons
  ps <- survminer::pairwise_survdiff(
    survival::Surv(time, status) ~ group3,
    data = input_pdata,
    p.adjust.method = "none"
  )

  addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001,
    "<0.001", round(ps$p.value, 3)
  )))
  addTab[is.na(addTab)] <- "-"

  df <- tibble::tibble(x = 0, y = 0, tb = list(addTab))

  rlang::check_installed("gridExtra")
  tb_grob <- gridExtra::tableGrob(df$tb,
    rows = TRUE,
    theme = gridExtra::ttheme_minimal(base_size = 6)
  )

  pp2$plot <- pp2$plot +
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$y, label = ""),
      data = df, size = 0
    ) +
    ggplot2::annotation_custom(tb_grob,
      xmin = df$x, xmax = df$x,
      ymin = df$y, ymax = df$y
    )

  res2 <- survminer::arrange_ggsurvplots(list(pp2), print = FALSE, ncol = 1, nrow = 1)

  if (!is.null(save_path)) {
    ggplot2::ggsave(
      plot = res2,
      filename = paste0(
        index, "-2-KMplot-3group-", signature, "-",
        project, ".", fig.type
      ),
      width = 6, height = 6.5, path = save_path
    )
  }

  # Plot 3: Two groups (median)
  input_pdata <- input_pdata[!is.na(input_pdata$group2), ]
  input_pdata$group2 <- ifelse(input_pdata$group2 == "High", 1, 0)

  cox_fit <- survival::coxph(
    survival::Surv(input_pdata$time, input_pdata$status) ~ group2,
    data = input_pdata
  )
  pvalue <- getHRandCIfromCoxph(cox_fit)
  HR <- paste0("Hazard Ratio = ", round(pvalue[, 2], 2))
  CI <- paste0("95% CI: ", round(pvalue[, 3], 2), " - ", round(pvalue[, 4], 2))

  sfit <- survminer::surv_fit(
    survival::Surv(input_pdata$time, input_pdata$status) ~ group2,
    data = input_pdata
  )

  pp3 <- survminer::ggsurvplot(
    sfit,
    data = input_pdata, censor = TRUE, ncensor.plot = FALSE,
    conf.int = FALSE, xlim = c(0, max_month), break.time.by = break_month_val,
    xlab = "Months after diagnosis",
    legend.labs = c(paste0("Low ", mini_sig), paste0("High ", mini_sig)),
    submain = paste0(signature, "-in-", project),
    surv.median.line = "h", risk.table = TRUE, tables.height = 0.20,
    palette = cols, pval.size = 8,
    pval = paste(
      ifelse(pvalue[, 1] < 0.0001, "P < 0.0001",
        paste0("P = ", round(pvalue[, 1], 4))
      ),
      HR, CI,
      sep = "\n"
    ),
    size = 0.4
  )

  res3 <- survminer::arrange_ggsurvplots(list(pp3), print = FALSE, ncol = 1, nrow = 1)

  if (!is.null(save_path)) {
    ggplot2::ggsave(
      plot = res3,
      filename = paste0(
        index, "-3-KMplot-2group-", signature, "-",
        project, ".", fig.type
      ),
      width = 6, height = 6.5, path = save_path
    )
  }

  # Combine all plots
  input_pdata$group2 <- ifelse(input_pdata$group2 == 1, "High", "Low")

  plots <- list(pp1, pp3, pp2)
  combined_plots <- survminer::arrange_ggsurvplots(plots, print = FALSE, ncol = 3, nrow = 1)

  list(data = input_pdata, plots = combined_plots)
}

#' Calculate break points for survival plots
#' @keywords internal
#' @noRd
calculate_break_month <- function(input, block = 6) {
  max_follow_up <- max(input, na.rm = TRUE)
  max_follow_up <- ceiling(max_follow_up)
  break_month_val <- max_follow_up / block
  break_month_val <- ceiling(break_month_val)
  break_month_val
}
