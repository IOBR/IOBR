#' Batch Signature Survival Plot
#'
#' @description
#' Generates Kaplan-Meier survival plots for multiple projects or cohorts based
#' on signature scores. Automatically determines optimal cutpoints for signature
#' stratification and creates publication-ready survival curves.
#'
#' @param input_pdata Data frame containing survival data and signature scores.
#' @param signature Character string specifying the column name of the target
#'   signature for survival analysis.
#' @param id Character string specifying the column name containing unique
#'   identifiers. Default is `"ID"`.
#' @param column_of_project Character string specifying the column name
#'   containing project identifiers. Default is `"ProjectID"`.
#' @param project Character string or vector specifying project name(s) to
#'   analyze. If `NULL`, all projects are analyzed. Default is `NULL`.
#' @param time Character string specifying the column name containing
#'   time-to-event data. Default is `"time"`.
#' @param status Character string specifying the column name containing event
#'   status. Default is `"status"`.
#' @param time_type Character string specifying the time unit. Options are
#'   `"day"` or `"month"`. Default is `"day"`.
#' @param break_month Numeric value or `"auto"` specifying the interval for
#'   time axis breaks in months. Default is `"auto"`.
#' @param palette Character string specifying the color palette. Default is
#'   `"jama"`.
#' @param cols Character vector of custom colors. If `NULL`, palette is used.
#'   Default is `NULL`.
#' @param mini_sig Character string for the signature label in the legend.
#'   Default is `"score"`.
#' @param save_path Character string specifying the directory path for saving
#'   plots.
#' @param show_col Logical indicating whether to display color information.
#'   Default is `TRUE`.
#' @param fig_type Character string specifying the output file format
#'   (`"pdf"`, `"png"`, etc.). Default is `"pdf"`.
#'
#' @return Data frame containing combined survival analysis results from all
#'   projects.
#'
#' @author Dongqiang Zeng
#' @export
#' @importFrom survminer surv_cutpoint
#' @import survival
#'
#' @examples
#' \dontrun{
#' sig_stad <- load_data("sig_stad")
#' result <- batch_sig_surv_plot(
#'   input_pdata = sig_stad,
#'   signature = "T.cells.CD8",
#'   id = "ID",
#'   column_of_project = "ProjectID",
#'   project = NULL,
#'   time = "OS_time",
#'   status = "OS_status",
#'   time_type = "month",
#'   break_month = "auto",
#'   palette = "jama",
#'   cols = NULL,
#'   mini_sig = "score",
#'   show_col = TRUE,
#'   fig_type = "pdf"
#' )
#' }
batch_sig_surv_plot <- function(
  input_pdata,
  signature,
  id = "ID",
  column_of_project = "ProjectID",
  project = NULL,
  time = "time",
  status = "status",
  time_type = "day",
  break_month = "auto",
  palette = "jama",
  cols = NULL,
  mini_sig = "score",
  save_path = file.path(tempdir(), "Multiple-KM-plot"),
  show_col = TRUE,
  fig_type = "pdf"
) {
  if (!is.data.frame(input_pdata)) {
    cli::cli_abort("{.arg input_pdata} must be a data frame")
  }
  if (!signature %in% colnames(input_pdata)) {
    cli::cli_abort("Signature column {.val {signature}} not found in input_pdata")
  }
  if (!column_of_project %in% colnames(input_pdata)) {
    cli::cli_abort("Project column {.val {column_of_project}} not found in input_pdata")
  }

  input_pdata <- as.data.frame(input_pdata)
  colnames(input_pdata)[colnames(input_pdata) == column_of_project] <- "ProjectID"
  input_pdata$ProjectID <- as.character(input_pdata$ProjectID)

  goi <- if (is.null(project)) {
    unique(input_pdata$ProjectID)
  } else {
    intersect(project, unique(input_pdata$ProjectID))
  }

  if (length(goi) == 0) {
    cli::cli_abort("No valid projects found")
  }

  results <- lapply(seq_along(goi), function(i) {
    var <- goi[i]
    cli::cli_alert_info("Processing project: {.val {var}}")

    pd <- input_pdata[input_pdata$ProjectID == var, , drop = FALSE]
    if (nrow(pd) == 0) {
      cli::cli_alert_warning("No data for project {.val {var}}, skipping")
      return(NULL)
    }

    sig_surv_plot(
      input_pdata = pd,
      signature = signature,
      project = var,
      ID = id,
      time = time,
      status = status,
      time_type = time_type,
      break_month = break_month,
      cols = cols,
      palette = palette,
      show_col = show_col,
      mini_sig = mini_sig,
      fig.type = fig_type,
      save_path = save_path,
      index = i
    )$data
  })

  invisible(dplyr::bind_rows(results))
}
