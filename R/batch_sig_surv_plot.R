#' Batch Signature Survival Plot
#'
#' @description
#' Generates Kaplan-Meier survival plots for multiple projects or cohorts based on
#' signature scores. Automatically determines optimal cutpoints for signature
#' stratification and creates publication-ready survival curves.
#'
#' @param input_pdata Data frame containing survival data and signature scores.
#' @param signature Character string specifying the column name of the target signature
#'   for survival analysis.
#' @param ID Character string specifying the column name containing unique identifiers.
#'   Default is \code{"ID"}.
#' @param column_of_Project Character string specifying the column name containing
#'   project identifiers. Default is \code{"ProjectID"}.
#' @param project Character string or vector specifying project name(s) to analyze.
#'   If \code{NULL}, all projects are analyzed. Default is \code{NULL}.
#' @param time Character string specifying the column name containing time-to-event data.
#'   Default is \code{"time"}.
#' @param status Character string specifying the column name containing event status.
#'   Default is \code{"status"}.
#' @param time_type Character string specifying the time unit. Options are \code{"day"}
#'   or \code{"month"}. Default is \code{"day"}.
#' @param break_month Numeric value or \code{"auto"} specifying the interval for time
#'   axis breaks in months. Default is \code{"auto"}.
#' @param palette Character string specifying the color palette. Default is \code{"jama"}.
#' @param cols Character vector of custom colors. If \code{NULL}, palette is used.
#'   Default is \code{NULL}.
#' @param mini_sig Character string for the signature label in the legend. Default is
#'   \code{"score"}.
#' @param save_path Character string specifying the directory path for saving plots.
#'   Default is \code{"Multiple-KM-plot"}.
#' @param show_col Logical indicating whether to display color information. Default is
#'   \code{TRUE}.
#' @param fig.type Character string specifying the output file format (\code{"pdf"},
#'   \code{"png"}, etc.). Default is \code{"pdf"}.
#'
#' @return Data frame containing combined survival analysis results from all projects.
#'
#' @author Dongqiang Zeng
#' @export
#' @importFrom survminer surv_cutpoint
#' @import survival
#' @examples
#' # Load TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' sig_stad <- as.data.frame(sig_stad)
#' # Generate survival plots for multiple projects
#' result <- batch_sig_surv_plot(
#'   input_pdata = sig_stad, signature = "T.cells.CD8", ID = "ID",
#'   column_of_Project = "ProjectID", project = NULL, time = "OS_time",
#'   status = "OS_status", time_type = "month", break_month = "auto",
#'   palette = "jama", cols = NULL, mini_sig = "score",
#'   save_path = "Multiple-KM-plot", show_col = TRUE, fig.type = "pdf"
#' )
batch_sig_surv_plot <- function(input_pdata, signature, ID = "ID", column_of_Project = "ProjectID", project = NULL, time = "time", status = "status", time_type = "day", break_month = "auto",
                                palette = "jama", cols = NULL, mini_sig = "score", save_path = "Multiple-KM-plot", show_col = TRUE, fig.type = "pdf") {
  input_pdata <- as.data.frame(input_pdata)
  goi <- as.character(levels(as.factor(as.character(input_pdata[, column_of_Project]))))

  # message(summary(as.factor(as.character(input_pdata[,column_of_Project]))))

  colnames(input_pdata)[which(colnames(input_pdata) == column_of_Project)] <- "ProjectID"

  input_pdata[, "ProjectID"] <- as.character(input_pdata[, "ProjectID"])
  ####################################

  input_pdata_com <- data.frame(NULL)

  for (i in 1:length(goi)) {
    var <- goi[i]
    message(">>> Preprocessing dataset: ", var)
    pd <- input_pdata[input_pdata$ProjectID %in% var, ]

    pd <- sig_surv_plot(
      input_pdata = pd,
      signature = signature,
      project = var,
      ID = ID,
      time = time,
      status = status,
      time_type = time_type,
      break_month = break_month,
      cols = cols,
      palette = palette,
      show_col = show_col,
      mini_sig = mini_sig,
      fig.type = fig.type,
      save_path = save_path,
      index = i
    )

    input_pdata_com <- rbind(input_pdata_com, pd$data)
  }

  return(input_pdata_com)
}
