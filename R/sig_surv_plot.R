#' Generate Kaplan-Meier Survival Plot for Signature
#'
#' Creates Kaplan-Meier survival plots for a given signature or gene, with automatic cutoff determination.
#'
#' @param input_pdata Data frame with survival data and signature scores.
#' @param signature Column name of the target signature.
#' @param project Project name for output. Default is "KM".
#' @param time Column name for survival time. Default is "time".
#' @param status Column name for survival status. Default is "status".
#' @param time_type Time unit. Default is "month".
#' @param break_month Time axis breaks. Default is "auto".
#' @param palette Color palette if cols not provided. Default is "jama".
#' @param save_path Directory for saving plots. Default is "KM-plot".
#' @param mini_sig Label for low score group. Default is "score".
#' @param show_col Logical indicating whether to show colors. Default is TRUE.
#' @param index Index for multiple plots. Default is 1.
#' @param cols Optional color vector.
#' @param fig.type File format. Default is "png".
#' @param ID Column name for sample IDs. Default is "ID".
#'
#' @return Survival analysis results.
#' @importFrom survminer surv_cutpoint
#' @import survival
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("tcga_stad_pdata", package = "IOBR")
#' sig_surv_plot(input_pdata = tcga_stad_pdata, signature = "TMEscore_plus",
#'               time = "time", status = "OS_status")
#'
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
                          save_path = "KM-plot",
                          index = 1) {
  if (!signature %in% colnames(input_pdata)) stop("Target signature must be in the column name of input_pdata")
  ###############################

  if (!file.exists(save_path)) dir.create(save_path)
  abspath <- paste(getwd(), "/", save_path, "/", sep = "")
  ###########################################

  # if("time"%in%input_pdata & time!="time") stop(">>>- 'time' is already in column name of input_pdata... ")
  # if("status"%in%input_pdata & time!="status") stop(">>>- 'status' is already in column name of input_pdata... ")

  input_pdata <- as.data.frame(input_pdata[, colnames(input_pdata) %in% c(ID, time, status, signature)])

  # filter data
  colnames(input_pdata)[which(colnames(input_pdata) == time)] <- "time"
  colnames(input_pdata)[which(colnames(input_pdata) == status)] <- "status"
  colnames(input_pdata)[which(colnames(input_pdata) == ID)] <- "ID"

  input_pdata$time <- as.numeric(input_pdata$time)
  input_pdata$status <- as.numeric(as.character(input_pdata$status))
  input_pdata <- input_pdata %>%
    filter(!is.na(time)) %>%
    filter(!is.na(status)) %>%
    filter(!.$time == "NA") %>%
    filter(!.$status == "NA") %>%
    filter(.$time > 0)
  ###################################
  input_pdata$time <- as.numeric(input_pdata$time)
  input_pdata$status <- as.numeric(input_pdata$status)
  ###################################

  input_pdata$ID <- as.character(input_pdata$ID)
  # input_pdata$ProjectID<-project
  ##################################

  # transform follow up time
  if (time_type == "day") input_pdata$time <- input_pdata$time / 30
  ###################################
  # cut use quantile
  #################################
  # print(paste0(signature,"'s expression range: "))
  # print( summary(cut(input_pdata[,signature],breaks = 3)))
  q1 <- quantile(input_pdata[, signature], probs = 1 / 3)
  q2 <- quantile(input_pdata[, signature], probs = 2 / 3)

  input_pdata$group3 <- input_pdata[, signature]
  input_pdata$group3 <- ifelse(input_pdata$group3 <= q1, "Low", ifelse(input_pdata$group3 >= q2, "High", "Middle"))

  input_pdata$group2 <- input_pdata[, signature]
  input_pdata$group2 <- ifelse(input_pdata$group2 <= mean(input_pdata[, signature]), "Low", "High")
  ####################################

  message(paste(
    ">>> Dataset's survival follow up time is range between",
    paste(round(summary(input_pdata$time), 2)[c(1, 6)], collapse = " to "), "months"
  ))
  ######################################
  # calculate best cutoff
  res.cut <- surv_cutpoint(input_pdata, time = "time", event = "status", variables = signature)
  res.cut <- res.cut$cutpoint[[1]]
  message(paste0(">>> The best cutoff for ", signature, " is: ", round(res.cut, 2)))
  ###################################
  # data category
  input_pdata <- best_cutoff(pdata = input_pdata, variable = signature, time = "time", status = "status", PrintResult = F)

  # input_pdata<-as.data.frame(input_pdata)
  #####################################
  # signatureBinary <- as.character(paste(signature,"_binary",sep = ""))
  colnames(input_pdata)[which(colnames(input_pdata) == paste0(signature, "_binary"))] <- "bestcutoff"


  message(paste(">>> High ", signature, " = ", summary(as.factor(input_pdata$bestcutoff)))[1], collapse = " ")
  message(paste(">>> Low ", signature, " = ", summary(as.factor(input_pdata$bestcutoff)))[2], collapse = " ")

  save(input_pdata, file = paste0(abspath, index, "-0-", project, "-", signature, "-survival-analysis-input.RData"))
  #######################################

  # Sur <-   Surv(time = input_pdata$time,event = input_pdata$status)
  # 将high变量转化为1，才能准确计算出HR，否则是反的
  #####################################
  input_pdata$bestcutoff <- ifelse(input_pdata$bestcutoff == "High", 1, 0)

  pvalue <- getHRandCIfromCoxph(coxph(Surv(time = input_pdata$time, event = input_pdata$status) ~ input_pdata$bestcutoff, data = input_pdata))

  HR <- paste("Hazard Ratio = ", round(pvalue[, 2], 2), sep = "")
  CI <- paste("95% CI: ", paste(round(pvalue[, 3], 2), round(pvalue[, 4], 2), sep = " - "), sep = "")
  cut_off <- paste("cutoff = ", round(res.cut, 3), sep = "")
  ###########################################
  sfit <- surv_fit(Surv(time = input_pdata$time, event = input_pdata$status) ~ input_pdata$bestcutoff,
    data = input_pdata
  )
  # input_pdata[,index]<-ifelse(input_pdata[,index]==1,"High","LOW")
  ############################################
  # define break time

  if (break_month == "auto") {
    break_month <- break_month(input = input_pdata$time, block = 6)
  }

  max_month <- break_month * 6
  ###########################################

  if (is.null(cols)) cols <- palettes(category = "box", palette = palette, show_col = F)

  ###########################################
  pp1 <- ggsurvplot(sfit,
    data = input_pdata,
    censor = TRUE,
    ncensor.plot = F, conf.int = F,
    xlim = c(0, max_month),
    break.time.by = break_month,
    xlab = "Months after diagnosis",
    surv.median.line = "h", # draw horizontal line at median survival
    legend.labs = c(paste0("Low ", mini_sig), paste0("High ", mini_sig)),
    submain = paste0(signature, "-in-", project),
    risk.table = T, tables.height = 0.20,
    palette = cols,
    pval.size = 8,
    pval = paste(
      pval = ifelse(pvalue[, 1] < 0.0001, "P < 0.0001",
        paste("P = ", round(pvalue[, 1], 4), sep = "")
      ),
      HR, CI, cut_off, sep = "\n"
    ),
    size = 0.4
  )

  res1 <- arrange_ggsurvplots(list(pp1), print = FALSE, ncol = 1, nrow = 1)
  ##############################

  ggsave(
    plot = res1, filename = paste0(index, "-1-KMplot-best-cutoff-", signature, "-", project, ".", fig.type),
    width = 6, height = 6.5, path = save_path
  )

  ################################################

  input_pdata$bestcutoff <- ifelse(input_pdata$bestcutoff == 1, "High", "Low")
  print(head(input_pdata))
  ################################################

  # choose  tri-sectional quantiles to conduct KM-plot
  ###########################################
  # turn high to 1 to calculate HR

  input_pdata <- input_pdata[!is.na(input_pdata$group3), ]
  # input_pdata$group3<-ifelse(input_pdata$group3=="High",3,ifelse(input_pdata$group3=="Middle",2,1))
  # Sur <- Surv(input_pdata$time,input_pdata$status)
  # pvalue<-getHRandCIfromCoxph(coxph(Surv(input_pdata$time,input_pdata$status)~input_pdata$group3,data = input_pdata))
  # HR <- paste("Hazard Ratio = ", round(pvalue[,2],2), sep = "")
  # CI <- paste("95% CI: ", paste(round(pvalue[,3],2), round(pvalue[,4],2), sep = " - "), sep = "")
  ##########################################
  sfit <- surv_fit(Surv(input_pdata$time, input_pdata$status) ~ input_pdata$group3, data = input_pdata)

  # hack strata for better survival curve
  names(sfit$strata) <- gsub("group3=", "", names(sfit$strata))
  ##########################################
  # ###########################################
  pp2 <- ggsurvplot(sfit,
    data = input_pdata, censor = TRUE,
    ncensor.plot = F, conf.int = F,
    xlim = c(0, max_month),
    break.time.by = break_month,
    xlab = "Months after diagnosis",
    # legend.labs = c(paste0('Low ',mini_sig),paste0("Middle ",mini_sig),paste0("High ",mini_sig)),
    submain = paste0(signature, "-in-", project),
    surv.median.line = "h", # draw horizontal line at median survival
    risk.table = T, tables.height = 0.25,
    palette = cols,
    pval.size = 8,
    size = 0.4
  )


  fitd <- survdiff(Surv(time, status) ~ group3,
    data = input_pdata,
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

  pp2$plot <- pp2$plot + annotate("text",
    x = 0, y = 0.55,
    hjust = 0,
    fontface = 3,
    label = p.lab
  )

  # calculate pair-wise survival comparison
  ps <- pairwise_survdiff(Surv(time, status) ~ group3,
    data = input_pdata,
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
  pp2$plot <- pp2$plot + ggpp::geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)


  res2 <- arrange_ggsurvplots(list(pp2), print = FALSE, ncol = 1, nrow = 1)
  ##############################

  ggsave(
    plot = res2, filename = paste0(index, "-2-KMplot-3group-", signature, "-", project, ".", fig.type),
    width = 6, height = 6.5, path = save_path
  )

  ################################################
  # input_pdata$group3<-ifelse(input_pdata$group3==3,"High",ifelse(input_pdata$group3==2,"Middle","Low"))
  ################################################

  # choose mean value-category data to conduct KM-plot
  # input_pdata = input_pdata
  # turn high to 1 to calculate HR
  input_pdata <- input_pdata[!is.na(input_pdata$group2), ]
  input_pdata$group2 <- ifelse(input_pdata$group2 == "High", 1, 0)
  # ###########################################
  # # turn high to 1 to calculate HR

  pvalue <- getHRandCIfromCoxph(coxph(Surv(input_pdata$time, input_pdata$status) ~ input_pdata[, "group2"], data = input_pdata))
  HR <- paste("Hazard Ratio = ", round(pvalue[, 2], 2), sep = "")
  CI <- paste("95% CI: ", paste(round(pvalue[, 3], 2), round(pvalue[, 4], 2), sep = " - "), sep = "")

  sfit <- surv_fit(Surv(input_pdata$time, input_pdata$status) ~ input_pdata$group2, data = input_pdata)
  pp3 <- ggsurvplot(sfit,
    data = input_pdata, censor = TRUE,
    ncensor.plot = F, conf.int = F,
    xlim = c(0, max_month),
    break.time.by = break_month,
    xlab = "Months after diagnosis",
    legend.labs = c(paste0("Low ", mini_sig), paste0("High ", mini_sig)),
    submain = paste0(signature, "-in-", project),
    surv.median.line = "h", # draw horizontal line at median survival
    risk.table = T,
    tables.height = 0.20,
    palette = cols,
    pval.size = 8,
    pval = paste(
      pval = ifelse(pvalue[, 1] < 0.0001, "P < 0.0001",
        paste("P = ", round(pvalue[, 1], 4), sep = "")
      ),
      HR, CI, sep = "\n"
    ),
    size = 0.4
  )

  res3 <- arrange_ggsurvplots(list(pp3), print = FALSE, ncol = 1, nrow = 1)
  ##############################
  ggsave(
    plot = res3, filename = paste0(index, "-3-KMplot-2group-", signature, "-", project, ".", fig.type),
    width = 6, height = 6.5, path = save_path
  )
  #############################################
  input_pdata$group2 <- ifelse(input_pdata$group2 == 1, "High", "Low")

  plots <- list()
  plots[[1]] <- pp1
  plots[[2]] <- pp3
  plots[[3]] <- pp2
  print(">>>>>>>>>")
  plots <- arrange_ggsurvplots(plots, print = FALSE, ncol = 3, nrow = 1)

  res <- list("data" = input_pdata, "plots" = plots)
  return(res)
}
