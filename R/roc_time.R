#' Time-dependent ROC Curve for Survival Analysis
#'
#' @description
#' Generates a time-dependent Receiver Operating Characteristic (ROC) plot to evaluate
#' the predictive performance of one or more variables in survival analysis. Calculates
#' the Area Under the Curve (AUC) for each specified time point and variable, and
#' creates a multi-line ROC plot with annotated AUC values.
#'
#' @param input Data frame containing variables for analysis.
#' @param vars Character vector. Variable(s) to be evaluated.
#' @param time Character string. Name of the time variable. Default is `"time"`.
#' @param status Character string. Name of the status variable. Default is `"status"`.
#' @param time_point Integer or vector. Time point(s) for ROC analysis. Default is `12`.
#' @param time_type Character string. Time unit (`"day"` or `"month"`). Default is `"month"`.
#' @param palette Character string. Color palette for the plot. Default is `"jama"`.
#' @param cols Character vector or string. Color specification: `"normal"`, `"random"`,
#'   or custom color vector.
#' @param seed Integer. Random seed for reproducibility. Default is `1234`.
#' @param show_col Logical indicating whether to display the color palette.
#'   Default is `FALSE`.
#' @param path Character string or `NULL`. Path to save the plot. Default is `NULL`.
#' @param main Character string. Main title of the plot. Default is `"PFS"`.
#' @param index Integer. Index for plot filename. Default is `1`.
#' @param fig.type Character string. Output file type (e.g., `"pdf"`, `"png"`).
#'   Default is `"pdf"`.
#' @param width Numeric. Width of the plot. Default is `5`.
#' @param height Numeric. Height of the plot. Default is `5.2`.
#'
#' @return A ggplot object representing the time-dependent ROC plot.
#'
#' @export
#' @author Dongqiang Zeng
#' @import survival
#'
#' @examples
#' tcga_stad_sig <- load_data("tcga_stad_sig")
#' pdata_stad <- load_data("pdata_stad")
#' input <- merge(pdata_stad, tcga_stad_sig, by = "ID")
#' roc_time(
#'   input = input, vars = c("Pan_F_TBRs", "CD_8_T_effector", "Immune_Checkpoint"),
#'   time = "time", status = "OS_status", time_point = 12, path = NULL, main = "OS"
#' )
roc_time <- function(input, vars, time = "time", status = "status", time_point = 12,
                     time_type = "month", palette = "jama", cols = "normal",
                     seed = 1234, show_col = FALSE, path = NULL, main = "PFS",
                     index = 1, fig.type = "pdf", width = 5, height = 5.2) {
  rlang::check_installed("timeROC")
  rlang::check_installed("survival")

  # Ensure survival is attached for timeROC internal calls
  if (!"package:survival" %in% search()) {
    attachNamespace("survival")
    on.exit(detach("package:survival", character.only = TRUE, unload = FALSE), add = TRUE)
  }

  save_plot <- !is.null(path)

  if (save_plot) {
    folder_info <- creat_folder(path)
    save_path <- folder_info$folder_name
  }

  # Get colors
  if (is.null(cols)) cols <- "normal"
  cols <- get_cols(cols = cols, palette = palette, show_col = show_col, seed = seed)

  # Prepare data
  input <- as.data.frame(input)
  input <- input[, colnames(input) %in% c(time, status, vars)]
  colnames(input)[colnames(input) == time] <- "time"
  colnames(input)[colnames(input) == status] <- "status"

  # Validate and convert
  time_type <- tolower(time_type)
  unit_label <- match.arg(time_type, c("day", "month"))

  input$time <- as.numeric(input$time)
  input$status <- as.numeric(as.character(input$status))

  input <- input %>%
    dplyr::filter(!is.na(.data$time)) %>%
    dplyr::filter(!is.na(.data$status)) %>%
    dplyr::filter(.data$time > 0)

  if (nrow(input) == 0) {
    cli::cli_abort("No valid data after filtering")
  }

  cli::cli_alert_info("Time range: {paste(round(range(input$time, na.rm = TRUE), 2), collapse = ' to ')}")

  time_vec <- input$time
  status_vec <- input$status

  # Route to appropriate handler
  n_vars <- length(vars)
  if (n_vars == 1) {
    p <- .roc_single_var(input, vars, time_vec, status_vec, time_point, unit_label, cols, main)
  } else if (n_vars == 2) {
    p <- .roc_two_vars(input, vars, time_vec, status_vec, time_point, unit_label, cols, main)
  } else if (n_vars == 3) {
    p <- .roc_three_vars(input, vars, time_vec, status_vec, time_point, unit_label, cols, main)
  } else {
    cli::cli_abort("Only 1-3 variables are supported")
  }

  p <- p + design_mytheme(axis_angle = 0, hjust = 0.5, axis_title_size = 1.7)

  if (save_plot) {
    ggplot2::ggsave(p,
      filename = paste0(index, "-", main, "-ROC-time.", fig.type),
      width = width, height = height, path = save_path
    )
  }

  p
}

#' ROC for single variable with 3 time points
#' @keywords internal
#' @noRd
.roc_single_var <- function(input, vars, time, status, time_point, unit_label, cols, main) {
  event_time <- time[status == 1]

  if (length(time_point) == 1) {
    if (length(event_time) >= 5) {
      time_point <- round(as.numeric(stats::quantile(event_time, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)), 0)
    } else {
      cli::cli_warn("Too few events; using follow-up time quantiles")
      time_point <- round(as.numeric(stats::quantile(time, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)), 0)
    }
  }

  if (length(time_point) != 3) {
    cli::cli_abort("For single variable, time_point must resolve to length 3")
  }

  time_point <- as.numeric(time_point)
  max_follow_up <- max(time, na.rm = TRUE)

  if (any(time_point > max_follow_up)) {
    cli::cli_abort("time_point exceeds maximum follow-up time ({max_follow_up})")
  }

  # Calculate ROCs
  rocs <- lapply(time_point, function(tp) {
    timeROC::timeROC(
      T = time, delta = status, marker = input[[vars]],
      cause = 1, weighting = "marginal", times = tp, iid = TRUE
    )
  })

  aucs <- vapply(rocs, function(r) round(as.numeric(r$AUC)[length(r$AUC)], 2), numeric(1))

  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[1]]$FP[, 2], y = rocs[[1]]$TP[, 2]), color = cols[1]) +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[2]]$FP[, 2], y = rocs[[2]]$TP[, 2]), color = cols[2]) +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[3]]$FP[, 2], y = rocs[[3]]$TP[, 2]), color = cols[3]) +
    ggplot2::geom_line(ggplot2::aes(x = c(0, 1), y = c(0, 1)), color = "grey", linetype = "dashed") +
    ggplot2::theme_light() +
    ggplot2::annotate("text",
      x = 0.7, y = 0.35,
      label = paste("AUC of", time_point[1], unit_label, "=", aucs[1]), color = cols[1]
    ) +
    ggplot2::annotate("text",
      x = 0.7, y = 0.25,
      label = paste("AUC of", time_point[2], unit_label, "=", aucs[2]), color = cols[2]
    ) +
    ggplot2::annotate("text",
      x = 0.7, y = 0.15,
      label = paste("AUC of", time_point[3], unit_label, "=", aucs[3]), color = cols[3]
    ) +
    ggplot2::scale_x_continuous(name = "False Positive Rate") +
    ggplot2::scale_y_continuous(name = "True Positive Rate") +
    ggplot2::ggtitle(paste0(vars, ", ", main, " = ", paste0(time_point, collapse = ", "), " ", unit_label))
}

#' ROC for two variables
#' @keywords internal
#' @noRd
.roc_two_vars <- function(input, vars, time, status, time_point, unit_label, cols, main) {
  if (length(time_point) != 1) {
    cli::cli_abort("For 2+ variables, time_point must be a single value")
  }

  max_follow_up <- max(time, na.rm = TRUE)
  if (time_point > max_follow_up) {
    cli::cli_abort("time_point ({time_point}) exceeds maximum follow-up ({max_follow_up})")
  }

  var1 <- vars[1]
  var2 <- vars[2]

  cox_fit <- survival::coxph(survival::Surv(time, status) ~ input[[var1]] + input[[var2]], data = input)
  combined_score <- stats::predict(cox_fit, type = "lp", newdata = input)

  markers <- list(input[[var1]], input[[var2]], combined_score)
  rocs <- lapply(markers, function(m) {
    timeROC::timeROC(
      T = time, delta = status, marker = m, cause = 1,
      weighting = "marginal", times = time_point, iid = TRUE
    )
  })

  aucs <- vapply(rocs, function(r) round(as.numeric(r$AUC)[length(r$AUC)], 2), numeric(1))

  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[1]]$FP[, 2], y = rocs[[1]]$TP[, 2]), color = cols[1]) +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[2]]$FP[, 2], y = rocs[[2]]$TP[, 2]), color = cols[2]) +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[3]]$FP[, 2], y = rocs[[3]]$TP[, 2]), color = cols[3]) +
    ggplot2::geom_line(ggplot2::aes(x = c(0, 1), y = c(0, 1)), color = "grey", linetype = "dashed") +
    ggplot2::theme_light() +
    ggplot2::annotate("text",
      x = 0.7, y = 0.35,
      label = paste("AUC of", var1, "=", aucs[1]), color = cols[1]
    ) +
    ggplot2::annotate("text",
      x = 0.7, y = 0.25,
      label = paste("AUC of", var2, "=", aucs[2]), color = cols[2]
    ) +
    ggplot2::annotate("text",
      x = 0.7, y = 0.15,
      label = paste("AUC of combined =", aucs[3]), color = cols[3]
    ) +
    ggplot2::scale_x_continuous(name = "False Positive Rate") +
    ggplot2::scale_y_continuous(name = "True Positive Rate") +
    ggplot2::ggtitle(paste0(main, " = ", time_point, " ", unit_label))
}

#' ROC for three variables
#' @keywords internal
#' @noRd
.roc_three_vars <- function(input, vars, time, status, time_point, unit_label, cols, main) {
  if (length(time_point) != 1) {
    cli::cli_abort("For 2+ variables, time_point must be a single value")
  }

  max_follow_up <- max(time, na.rm = TRUE)
  if (time_point > max_follow_up) {
    cli::cli_abort("time_point ({time_point}) exceeds maximum follow-up ({max_follow_up})")
  }

  var1 <- vars[1]
  var2 <- vars[2]
  var3 <- vars[3]

  cox_fit <- survival::coxph(survival::Surv(time, status) ~ input[[var1]] + input[[var2]] + input[[var3]], data = input)
  combined_score <- stats::predict(cox_fit, type = "lp", newdata = input)

  markers <- list(input[[var1]], input[[var2]], input[[var3]], combined_score)
  rocs <- lapply(markers, function(m) {
    timeROC::timeROC(
      T = time, delta = status, marker = m, cause = 1,
      weighting = "marginal", times = time_point, iid = TRUE
    )
  })

  aucs <- vapply(rocs, function(r) round(as.numeric(r$AUC)[length(r$AUC)], 2), numeric(1))

  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[1]]$FP[, 2], y = rocs[[1]]$TP[, 2]), color = cols[1]) +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[2]]$FP[, 2], y = rocs[[2]]$TP[, 2]), color = cols[2]) +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[3]]$FP[, 2], y = rocs[[3]]$TP[, 2]), color = cols[3]) +
    ggplot2::geom_line(ggplot2::aes(x = rocs[[4]]$FP[, 2], y = rocs[[4]]$TP[, 2]), color = cols[4]) +
    ggplot2::geom_line(ggplot2::aes(x = c(0, 1), y = c(0, 1)), color = "grey", linetype = "dashed") +
    ggplot2::theme_light() +
    ggplot2::annotate("text",
      x = 0.7, y = 0.35,
      label = paste("AUC of", var1, "=", aucs[1]), color = cols[1]
    ) +
    ggplot2::annotate("text",
      x = 0.7, y = 0.25,
      label = paste("AUC of", var2, "=", aucs[2]), color = cols[2]
    ) +
    ggplot2::annotate("text",
      x = 0.7, y = 0.15,
      label = paste("AUC of", var3, "=", aucs[3]), color = cols[3]
    ) +
    ggplot2::annotate("text",
      x = 0.7, y = 0.05,
      label = paste("AUC of combined =", aucs[4]), color = cols[4]
    ) +
    ggplot2::scale_x_continuous(name = "False Positive Rate") +
    ggplot2::scale_y_continuous(name = "True Positive Rate") +
    ggplot2::ggtitle(paste0(main, " = ", time_point, " ", unit_label))
}
