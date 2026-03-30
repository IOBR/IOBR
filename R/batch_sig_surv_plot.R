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
#'   plots. Default is `"Multiple-KM-plot"`.
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
#' # Load TCGA-STAD microenvironment signature data
#' sig_stad <- load_data("sig_stad")
#' sig_stad <- as.data.frame(sig_stad)
#' # Generate survival plots for multiple projects
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
#'   save_path = "Multiple-KM-plot",
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
    save_path = "Multiple-KM-plot",
    show_col = TRUE,
    fig_type = "pdf") {
  input_pdata <- as.data.frame(input_pdata)
  goi <- as.character(levels(as.factor(
    as.character(input_pdata[, column_of_project])
  )))

  colnames(input_pdata)[
    which(colnames(input_pdata) == column_of_project)
  ] <- "ProjectID"

  input_pdata[, "ProjectID"] <- as.character(input_pdata[, "ProjectID"])

  input_pdata_com <- data.frame(NULL)

  for (i in seq_along(goi)) {
    var <- goi[i]
    message(">>> Preprocessing dataset: ", var)
    pd <- input_pdata[input_pdata$ProjectID %in% var, ]

    pd <- sig_surv_plot(
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
    )

    input_pdata_com <- rbind(input_pdata_com, pd$data)
  }

  input_pdata_com
}
