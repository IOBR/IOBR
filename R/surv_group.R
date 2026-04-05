#' Generate Kaplan-Meier Survival Plots for Categorical Groups
#'
#' @description
#' Creates Kaplan-Meier survival plots for data grouped by a categorical variable.
#' Handles both binary and multi-level categorical groups with customizable plot
#' aesthetics.
#'
#' @param input_pdata Data frame containing survival data and grouping variables.
#' @param target_group Name of column containing the grouping variable.
#' @param ID Name of column with unique identifiers. Default is "ID".
#' @param levels Names for levels of target_group (for binary groups).
#'   Default is c("High", "Low").
#' @param reference_group Reference level for binary comparison. Default is NULL.
#' @param project Optional title for plot. Default is NULL.
#' @param time Name of column with follow-up times. Default is "time".
#' @param status Name of column with event indicators. Default is "status".
#' @param time_type Units: "month" or "day". Default is "month".
#' @param break_month X-axis break interval. If "auto", calculated automatically.
#'   Default is "auto".
#' @param cols Color vector for plot lines. Default is NULL.
#' @param palette Color palette name. Default is "jama".
#' @param mini_sig Prefix label for variables. Default is "score".
#' @param save_path Directory for saving plot. Default is NULL.
#' @param fig.type File format: "pdf" or "png". Default is "pdf".
#' @param index Identifier for file naming. Default is 1.
#' @param width Plot width. Default is 6.
#' @param height Plot height. Default is 6.5.
#' @param font.size.table Font size for risk table. Default is 3.
#'
#' @return Kaplan-Meier plot object.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' tcga_stad_pdata <- load_data("tcga_stad_pdata")
#' \dontrun{
#' surv_group(
#'   input_pdata = tcga_stad_pdata,
#'   target_group = "TMEscore_plus_binary",
#'   time = "time",
#'   status = "OS_status"
#' )
#' }
surv_group <- function(input_pdata,
                       target_group,
                       ID = "ID",
                       levels = c("High", "Low"),
                       reference_group = NULL,
                       project = NULL,
                       time = "time",
                       status = "status",
                       time_type = "month",
                       break_month = "auto",
                       cols = NULL,
                       palette = "jama",
                       mini_sig = "score",
                       save_path = NULL,
                       fig.type = "pdf",
                       index = 1,
                       width = 6,
                       height = 6.5,
                       font.size.table = 3) {
  if (!target_group %in% colnames(input_pdata)) {
    cli::cli_abort("{.arg target_group} must be a column in {.arg input_pdata}")
  }

  # Setup save path
  if (!is.null(save_path)) {
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    abspath <- file.path(normalizePath(save_path, winslash = "/", mustWork = FALSE), "")
  } else {
    abspath <- NULL
  }

  # Prepare data
  input_pd <- input_pdata[, c(ID, target_group, time, status)]
  colnames(input_pd)[which(colnames(input_pd) == time)] <- "time"
  colnames(input_pd)[which(colnames(input_pd) == status)] <- "status"
  colnames(input_pd)[which(colnames(input_pd) == ID)] <- "ID"

  # Convert and filter
  input_pd$time <- as.numeric(as.character(input_pd$time))
  input_pd$status <- as.numeric(as.character(input_pd$status))

  input_pd <- input_pd %>%
    dplyr::filter(!is.na(.data$time), !is.na(.data$status)) %>%
    dplyr::filter(.data$time > 0)

  # Transform time if needed
  if (time_type == "day") {
    input_pd$time <- input_pd$time / 30
  }

  cli::cli_alert_info(
    "Follow-up time ranges from {paste(round(summary(input_pd$time), 2)[c(1, 6)], collapse = ' to ')} months"
  )

  if (!is.null(save_path)) {
    save(input_pd, file = file.path(abspath, paste0("0-", project, "-", target_group, "-survival.RData")))
  }

  colnames(input_pd)[which(colnames(input_pd) == target_group)] <- "target_group"
  input_pd <- input_pd[!is.na(input_pd$target_group), ]
  print(summary(as.factor(input_pd$target_group)))

  # Define break time and colors
  break_month_val <- break_month
  if (identical(break_month_val, "auto")) {
    break_month_val <- calculate_break_month(input = input_pd$time, block = 6)
  }
  break_month_val <- as.numeric(break_month_val)
  max_month <- break_month_val * 6

  if (is.null(cols)) {
    cols <- palettes(category = "box", palette = palette, show_col = FALSE)
  }

  # Create plot based on number of groups
  n_groups <- length(unique(input_pd$target_group))

  if (n_groups > 2) {
    pp <- .surv_plot_multi(
      input_pd, max_month, break_month_val, cols,
      target_group, project, font.size.table
    )
  } else {
    pp <- .surv_plot_binary(
      input_pd, max_month, break_month_val, cols,
      target_group, project, reference_group, levels
    )
  }

  pp <- list(pp)
  res <- survminer::arrange_ggsurvplots(pp, print = FALSE, ncol = 1, nrow = 1)

  if (!is.null(save_path)) {
    ggplot2::ggsave(res,
      filename = paste0(index, "-KMplot-", target_group, "-", project, ".", fig.type),
      width = width, height = height, path = save_path
    )
  }

  res
}

#' @keywords internal
.surv_plot_multi <- function(input_pd, max_month, break_month, cols,
                             target_group, project, font.size.table) {
  sfit <- survminer::surv_fit(
    survival::Surv(time, status) ~ target_group,
    data = input_pd
  )

  names(sfit$strata) <- gsub("target_group=", "", names(sfit$strata))

  pp <- survminer::ggsurvplot(sfit,
    data = input_pd,
    censor = TRUE,
    ncensor.plot = FALSE,
    conf.int = FALSE,
    xlim = c(0, max_month),
    break.time.by = break_month,
    xlab = "Months after diagnosis",
    submain = paste0(target_group, " ", project),
    surv.median.line = "h",
    risk.table = TRUE,
    tables.height = 0.25,
    palette = cols,
    pval.size = 8
  )

  # Log-rank test
  fitd <- survival::survdiff(
    survival::Surv(time, status) ~ target_group,
    data = input_pd,
    na.action = stats::na.exclude
  )
  p.val <- 1 - stats::pchisq(fitd$chisq, length(fitd$n) - 1)

  p.lab <- paste0(
    "Overall P",
    ifelse(p.val < 0.001, " < 0.001", paste0(" = ", round(p.val, 3)))
  )

  pp$plot <- pp$plot + ggplot2::annotate("text",
    x = 0, y = 0.55,
    hjust = 0,
    fontface = 3,
    label = p.lab
  )

  # Pairwise comparison table
  ps <- survminer::pairwise_survdiff(
    survival::Surv(time, status) ~ target_group,
    data = input_pd,
    p.adjust.method = "none"
  )

  addTab <- as.data.frame(as.matrix(ifelse(
    round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3)
  )))
  addTab[is.na(addTab)] <- "-"

  df <- tibble::tibble(x = 0, y = 0, tb = list(addTab))
  rlang::check_installed("gridExtra")

  tb_grob <- gridExtra::tableGrob(df$tb,
    rows = TRUE,
    theme = gridExtra::ttheme_minimal(base_size = font.size.table)
  )

  pp$plot <- pp$plot +
    ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = ""), data = df, size = 0) +
    ggplot2::annotation_custom(tb_grob, xmin = df$x, xmax = df$x, ymin = df$y, ymax = df$y)

  pp
}

#' @keywords internal
.surv_plot_binary <- function(input_pd, max_month, break_month, cols,
                              target_group, project, reference_group, levels) {
  levels <- unique(input_pd$target_group)
  levels <- levels[order(levels)]

  if (!is.null(reference_group)) {
    if (!reference_group %in% input_pd$target_group) {
      cli::cli_abort("{.arg reference_group} must be one of target_group levels")
    }
    if (reference_group != levels[1]) {
      levels <- c(levels[2], levels[1])
    }
    input_pd$target_group <- ifelse(input_pd$target_group == reference_group, 1, 0)
  } else {
    cli::cli_alert_info("Reference group not defined, using alphabetical order")
  }

  pvalue <- getHRandCIfromCoxph(survival::coxph(
    survival::Surv(time, status) ~ target_group,
    data = input_pd
  ))

  HR <- paste("Hazard Ratio =", round(pvalue[, 2], 2))
  CI <- paste("95% CI:", paste(round(pvalue[, 3], 2), round(pvalue[, 4], 2), sep = " - "))

  sfit <- survminer::surv_fit(
    survival::Surv(time, status) ~ target_group,
    data = input_pd
  )

  survminer::ggsurvplot(sfit,
    data = input_pd,
    censor = TRUE,
    ncensor.plot = FALSE,
    conf.int = FALSE,
    xlim = c(0, max_month),
    break.time.by = break_month,
    xlab = "Months after diagnosis",
    legend.labs = c(levels[2], levels[1]),
    submain = paste0(target_group, " ", project),
    risk.table = TRUE,
    tables.height = 0.20,
    palette = cols,
    pval.size = 8,
    pval = paste(
      ifelse(pvalue[, 1] < 0.0001, "P < 0.0001", paste("P =", round(pvalue[, 1], 4))),
      HR, CI,
      sep = "\n"
    )
  )
}

#' Break Time Into Blocks
#'
#' @description
#' Divides time duration into specified blocks for analysis.
#'
#' @param input Numeric vector of time durations.
#' @param block Number of blocks. Default is 6.
#' @param time_type Units: "month" or "day". Default is "month".
#'
#' @return Numeric vector of breakpoints, rounded to nearest multiple of 5.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' time_data <- c(24, 36, 12, 48)
#' blocks <- calculate_break_month(input = time_data)
calculate_break_month <- function(input, block = 6, time_type = c("month", "day")) {
  time_type <- rlang::arg_match(time_type)

  max_time <- max(input, na.rm = TRUE)
  if (time_type == "day") {
    max_time <- max_time / 30
  }

  cli::cli_alert_info(
    "Maximum follow-up time is {round(max_time, 1)} months; divided into {block} sections"
  )

  round((max_time %/% block) / 5) * 5
}
