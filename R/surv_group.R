#' Generate Kaplan-Meier Survival Plots for Categorical Groups
#'
#' This function creates Kaplan-Meier survival plots for data grouped by a categorical variable.
#' It handles both binary and multi-level categorical groups, adjusts follow-up times, and allows for
#' customizable plot aesthetics and output formats.
#'
#' @param input_pdata A data frame containing survival data and grouping variables.
#' @param target_group The name of the column in input_pdata that contains the grouping variable.
#' @param ID The name of the column in input_pdata that serves as a unique identifier for each row.
#' @param levels A vector of names for the levels of the target_group; used for labeling in the plot.
#' @param reference_group The reference level for comparison in survival analysis; used only for binary groups.
#' @param project An optional title for the plot; used in file naming and plot subtitles.
#' @param time The name of the column containing follow-up times.
#' @param status The name of the column containing event indicators (1=event occurred, 0=censored).
#' @param time_type Units of follow-up time ('month' or 'day'). Converts days to months if 'day' is selected.
#' @param break_month The interval for breaking time on the X-axis of the plot. If 'auto', it is calculated automatically.
#' @param cols A vector of colors for the plot lines; if NULL, uses a default palette based on 'palette' parameter.
#' @param palette The name of the color palette to use if 'cols' is NULL. Common choices include 'jama' and other palette names from the `ggsci` package.
#' @param mini_sig A prefix label for variables in the plot, typically indicating a scoring metric.
#' @param save_path The directory path where the plot should be saved.
#' @param fig.type The file format for the saved plot ('pdf' or 'png').
#' @param index A unique identifier for the plot file, typically used in file naming.
#' @param width The width of the output plot.
#' @param height The height of the output plot.
#' @param font.size.table The font size used in the plot's risk table.
#'
#' @return Saves the Kaplan-Meier plot to the specified path and returns the plot object.
#' @import dplyr
#' @import ggplot2
#' @import survival
#' @import survminer
#' @author Dongqiang Zeng
#'
#' @examples
#' data("tcga_stad_pdata", package = "IOBR")
#' surv_group(input_pdata = tcga_stad_pdata, target_group = "TMEscore_plus_binary", time = "time", status = "OS_status")
#' @export
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
                       save_path = paste0("KMplot"),
                       fig.type = "pdf",
                       index = 1,
                       width = 6,
                       height = 6.5,
                       font.size.table = 3) {
  if (!target_group %in% colnames(input_pdata)) stop("Target group must be in the column of input_pdata")
  ###############################

  if (!file.exists(save_path)) dir.create(save_path)
  abspath <- paste(getwd(), "/", save_path, "/", sep = "")
  ###########################################
  input_pd <- input_pdata[, c(ID, target_group, time, status)]

  # if(time!="time" & "time"%in%colnames(input_pd)) stop(paste0("time already exists in column name of input_pd."))
  # if(status!="status" & "status"%in%colnames(input_pd)) stop(paste0("status already exists in column name of input_pd."))

  # filter data
  colnames(input_pd)[which(colnames(input_pd) == time)] <- "time"
  colnames(input_pd)[which(colnames(input_pd) == status)] <- "status"
  colnames(input_pd)[which(colnames(input_pd) == ID)] <- "ID"


  input_pd$time <- as.numeric(input_pd$time)
  input_pd$status <- as.numeric(as.character(input_pd$status))
  input_pd <- input_pd %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::filter(!is.na(status)) %>%
    dplyr::filter(!.$time == "NA") %>%
    dplyr::filter(!.$status == "NA") %>%
    dplyr::filter(.$time > 0)
  ###################################
  input_pd$time <- as.numeric(input_pd$time)
  input_pd$status <- as.numeric(input_pd$status)
  ###################################
  input_pd <- as.data.frame(input_pd[, colnames(input_pd) %in% c("ID", "time", "status", target_group)])
  input_pd$ID <- as.character(input_pd$ID)
  # input_pd$ProjectID<-project
  ##################################

  # transform follow up time
  if (time_type == "day") input_pd$time <- input_pd$time / 30
  ###################################
  # cut use quantile
  #################################
  message(paste(
    ">>> Dataset's survival follow up time is range between",
    paste(round(summary(input_pd$time), 2)[c(1, 6)], collapse = " to "), "months"
  ))


  ###################################
  save(input_pd, file = paste0(abspath, "0-", project, "-", target_group, "-survival-analysis-input.RData"))
  #######################################

  colnames(input_pd)[which(colnames(input_pd) == target_group)] <- "target_group"
  input_pd <- input_pd[!is.na(input_pd$target_group), ]
  message(print(summary(as.factor(input_pd$target_group))))
  ####################################

  # define the break time and colors
  if (break_month == "auto") {
    break_month <- break_month(input = input_pd$time, block = 6)
  }
  max_month <- break_month * 6
  ###########################################
  if (is.null(cols)) cols <- palettes(category = "box", palette = palette, show_col = F)


  # print(unique(input_pd$target_group))

  if (length(unique(input_pd$target_group)) > 2) {
    input_pd <- input_pd[!is.na(input_pd$target_group), ]
    # input_pd$target_group<-ifelse(input_pd$target_group=="High",3,ifelse(input_pd$target_group=="Middle",2,1))
    # Sur <- Surv(input_pd$time,input_pd$status)
    # pvalue<-getHRandCIfromCoxph(coxph(Surv(input_pd$time,input_pd$status)~input_pd$target_group,data = input_pd))
    # HR <- paste("Hazard Ratio = ", round(pvalue[,2],2), sep = "")
    # CI <- paste("95% CI: ", paste(round(pvalue[,3],2), round(pvalue[,4],2), sep = " - "), sep = "")
    ##########################################
    sfit <- survminer::surv_fit(Surv(input_pd$time, input_pd$status) ~ input_pd$target_group, data = input_pd)

    # hack strata for better survival curve
    names(sfit$strata) <- gsub("target_group=", "", names(sfit$strata))
    #############################################
    #############################################
    pp <- survminer::ggsurvplot(sfit,
      data = input_pd,
      censor = TRUE,
      ncensor.plot = F, conf.int = F,
      xlim = c(0, max_month),
      break.time.by = break_month,
      xlab = "Months after diagnosis",
      # legend.labs    = c(paste0('Low ',mini_sig),paste0("Middle ",mini_sig),paste0("High ",mini_sig)),
      submain = paste0(target_group, " ", project),
      surv.median.line = "h", # draw horizontal line at median survival
      risk.table = T,
      tables.height = 0.25,
      palette = cols,
      pval.size = 8
    )

    fitd <- survdiff(Surv(time, status) ~ target_group,
      data = input_pd,
      na.action = na.exclude
    )

    p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

    # add nominal pvalue for log-rank test
    p.lab <- paste0(
      "Overall P",
      ifelse(p.val < 0.001, " < 0.001",
        paste0(" = ", round(p.val, 3))
      )
    )

    pp$plot <- pp$plot + annotate("text",
      x = 0, y = 0.55,
      hjust = 0,
      fontface = 3,
      # size = 1,
      label = p.lab
    )

    # calculate pair-wise survival comparison
    ps <- pairwise_survdiff(Surv(time, status) ~ target_group,
      data = input_pd,
      p.adjust.method = "none"
    )

    # add pair-wise comparison table
    # options(stringsAsFactors = FALSE)
    addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
      round(ps$p.value, 3)
    )))
    addTab[is.na(addTab)] <- "-"
    # options(stringsAsFactors = TRUE)

    df <- tibble(x = 0, y = 0, tb = list(addTab))
    pp$plot <- pp$plot + ggpp::geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE, size = font.size.table)
  } else {
    if (is.null(reference_group)) {
      levels <- unique(input_pd$target_group)
      levels <- levels[order(levels)]
      #####################################
      message(paste0(">>>--- Reference_group was not defined..."))

      pvalue <- getHRandCIfromCoxph(coxph(Surv(time = input_pd$time, event = input_pd$status) ~ input_pd$target_group, data = input_pd))

      HR <- paste("Hazard Ratio = ", round(pvalue[, 2], 2), sep = "")
      CI <- paste("95% CI: ", paste(round(pvalue[, 3], 2), round(pvalue[, 4], 2), sep = " - "), sep = "")
      # cut_off<-paste("cutoff = ",round(res.cut,3),sep = "")
      ###########################################
      sfit <- surv_fit(Surv(time = input_pd$time, event = input_pd$status) ~ input_pd$target_group, data = input_pd)
      # input_pd[,index]<-ifelse(input_pd[,index]==1,"High","LOW")
      ###########################################
      pp <- survminer::ggsurvplot(sfit,
        data = input_pd,
        censor = TRUE,
        ncensor.plot = F, conf.int = F,
        xlim = c(0, max_month),
        break.time.by = break_month,
        xlab = "Months after diagnosis",
        legend.labs = c(paste0(levels[1]), paste0(levels[2])),
        submain = paste0(target_group, " ", project),
        risk.table = T,
        tables.height = 0.20,
        palette = cols,
        pval.size = 8,
        pval = paste(
          pval = ifelse(pvalue[, 1] < 0.0001, "P < 0.0001",
            paste("P = ", round(pvalue[, 1], 4), sep = "")
          ),
          HR, CI, sep = "\n"
        )
      )
    } else {
      # print(unique(input_pd$target_group))
      levels <- unique(input_pd$target_group)
      levels <- levels[order(levels)]
      if (reference_group != levels[1]) {
        levels <- c(levels[2], levels[1])
      }
      #####################################
      if (!reference_group %in% c(input_pd$target_group)) stop(">>>--- Reference_group must be one of target...")

      input_pd$target_group <- ifelse(input_pd$target_group == reference_group, 1, 0)

      pvalue <- getHRandCIfromCoxph(coxph(Surv(time = input_pd$time, event = input_pd$status) ~ input_pd$target_group, data = input_pd))

      HR <- paste("Hazard Ratio = ", round(pvalue[, 2], 2), sep = "")

      CI <- paste("95% CI: ", paste(round(pvalue[, 3], 2), round(pvalue[, 4], 2), sep = " - "), sep = "")
      # cut_off<-paste("cutoff = ",round(res.cut,3),sep = "")
      ###########################################
      sfit <- surv_fit(Surv(time = input_pd$time, event = input_pd$status) ~ input_pd$target_group, data = input_pd)
      # input_pd[,index]<-ifelse(input_pd[,index]==1,"High","LOW")
      ############################################
      # define break time
      ###########################################
      pp <- survminer::ggsurvplot(sfit,
        data = input_pd,
        censor = TRUE,
        ncensor.plot = F, conf.int = F,
        xlim = c(0, max_month),
        break.time.by = break_month,
        xlab = "Months after diagnosis",
        legend.labs = c(paste0(levels[2]), paste0(levels[1])),
        submain = paste0(target_group, " ", project),
        risk.table = T,
        tables.height = 0.20,
        palette = cols,
        pval.size = 8,
        pval = paste(
          pval = ifelse(pvalue[, 1] < 0.0001, "P < 0.0001",
            paste("P = ", round(pvalue[, 1], 4), sep = "")
          ),
          HR, CI, sep = "\n"
        )
      )
    }
  }
  pp <- list(pp)
  res <- arrange_ggsurvplots(pp, print = FALSE, ncol = 1, nrow = 1)

  ggsave(res,
    filename = paste0(index, "-KMplot-", target_group, "-", project, ".", fig.type),
    width = width, height = height, path = save_path
  )
  return(res)
}





#' Break Time Into Blocks
#'
#' Divides a time duration into specified blocks. This function is useful for creating intervals or categories
#' within a given time period, such as months or days, for further analysis or visualization in studies where
#' time segmentation might be relevant.
#'
#' @param input A numeric vector representing time durations that need to be divided, typically in months or days.
#' @param block An integer specifying the number of blocks to divide the time into; default is 6.
#' @param time_type A character string specifying the units of the input time: "month" for months and "day" for days.
#'        The default is "month". If "day" is specified, the function converts days into months by dividing by 30.
#'
#' @return A numeric vector representing the breakpoints for the time blocks, rounded to the nearest multiple of 5.
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' # Example with time in months
#' time_data <- c(24, 36, 12, 48)
#' blocks <- break_month(input = time_data)
#' print(blocks)
#'
#' # Example with time in days
#' day_data <- c(720, 1080, 360, 1440) # Corresponding to 24, 36, 12, 48 months
#' blocks_days <- break_month(input = day_data, time_type = "day")
#' print(blocks_days)
break_month <- function(input, block = 6, time_type = "month") {
  max_time <- max(input)
  if (time_type == "month") {
    max_time <- max_time
  } else if (time_type == "day") {
    max_time <- max_time / 30
  }

  message(paste0("  Maximum of follow up time is ", max_time, " months; and will be divided into ", block, " sections;"))
  round(c(max_time %/% block) / 5, 0) * 5
}
