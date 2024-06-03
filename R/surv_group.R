





#' Drawing KM-plot to a Categorical variable with two or more groups
#'
#' This function produces Kaplan-Meier survival plots for a given categorical variable with two or more groups. It allows customization of the analysis and visualization, including defining the time scale (months or days), breaking the follow-up time into intervals, selecting color palettes, and saving the plots. The function supports comparisons between groups using survival statistics and can handle additional grouping variables.
#'
#' @param input_pdata Data frame containing the survival time, status, and group variables.
#' @param target_group The name of the categorical variable in `input_pdata` used to define groups for survival analysis.
#' @param ID Identifier for the data rows; default is "ID".
#' @param levels A vector of the levels or groups within the target group variable; defaults to c("High", "Low").
#' @param reference_group The reference group against which other groups are compared in survival analysis; defaults to the first level if not specified.
#' @param project Name for the main title of the plot and part of the filename for saving the plot.
#' @param time The column name in `input_pdata` representing follow-up time; default is "time".
#' @param status The column name in `input_pdata` representing the event status; default is "status".
#' @param time_type Unit of follow-up time, either "month" or "day"; default is "month".
#' @param break_month Time intervals for breaking the follow-up period; default is "auto" which calculates based on data.
#' @param palette Color palette for plotting; default is "jama".
#' @param save_path Directory where the output files will be stored; default is "KMplot".
#' @param mini_sig A prefix used for variable naming in plots and files; default is "score".
#' @param cols Color vector for plots; if NULL, colors are selected based on `palette`.
#' @param fig.type The format for saving plots, either "pdf" or "png"; default is "pdf".
#' @param width Width of the output plot; default is 6 inches.
#' @param height Height of the output plot; default is 6.5 inches.
#' @param font.size.table Font size for tables included in the plot; default is 62.
#' @param index Folder name or prefix for saving the plot; used in file naming.
#'
#' @author Dongqiang Zeng
#' @import survival
#' @return A list containing the ggsurvplot object and optionally comparison results if comparisons are enabled.
#' @export
#'
#' @examples
#'data("tcga_stad_pdata", package = "IOBR")
#'surv_group(input_pdata = tcga_stad_pdata, target_group = "TMEscore_plus_binary", time = "time", status = "OS_status")
#'
surv_group  <-function(input_pdata,
                       target_group,
                       ID              = "ID",
                       levels          = c("High","Low"),
                       reference_group = "High",
                       project         = NULL,
                       time            = "time",
                       status          = "status",
                       time_type       = "month",
                       break_month     = "auto",
                       cols            = NULL,
                       palette         ="jama",
                       mini_sig        = "score",
                       save_path       = paste0("KMplot"),
                       fig.type        = "pdf",
                       index           = 1,
                       width           = 6,
                       height          = 6.5,
                       font.size.table = 3){


  if(!target_group%in%colnames(input_pdata)) stop("Target group must be in the column of input_pdata")
  ###############################

  if (!file.exists(save_path)) dir.create(save_path)
  abspath<-paste(getwd(),"/" ,save_path, "/",sep ="" )
  ###########################################
  input_pd <- input_pdata[,c(ID, target_group, time, status)]

  # if(time!="time" & "time"%in%colnames(input_pd)) stop(paste0("time already exists in column name of input_pd."))
  # if(status!="status" & "status"%in%colnames(input_pd)) stop(paste0("status already exists in column name of input_pd."))

  #filter data
  colnames(input_pd)[which(colnames(input_pd)==time)]<-"time"
  colnames(input_pd)[which(colnames(input_pd)==status)]<-"status"
  colnames(input_pd)[which(colnames(input_pd)==ID)]<-"ID"


  input_pd$time<-as.numeric(input_pd$time)
  input_pd$status<-as.numeric(as.character(input_pd$status))
  input_pd<- input_pd %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::filter(!is.na(status)) %>%
    dplyr::filter(!.$time=="NA") %>%
    dplyr::filter(!.$status=="NA") %>%
    dplyr::filter(.$time > 0)
  ###################################
  input_pd$time<-as.numeric(input_pd$time)
  input_pd$status<-as.numeric(input_pd$status)
  ###################################
  input_pd<-as.data.frame(input_pd[,colnames(input_pd)%in%c("ID","time","status",target_group)])
  input_pd$ID<-as.character(input_pd$ID)
  # input_pd$ProjectID<-project
  ##################################

  #transform follow up time
  if(time_type=="day") input_pd$time<-input_pd$time/30
  ###################################
  #cut use quantile
  #################################
  message(paste(">>> Dataset's survival follow up time is range between",
                paste(round(summary(input_pd$time),2)[c(1,6)],collapse = " to "),"months"))


  ###################################
  save(input_pd,file = paste0(abspath,"0-",project,"-",target_group,"-survival-analysis-input.RData"))
  #######################################

  colnames(input_pd)[which(colnames(input_pd)==target_group)]<-"target_group"
  input_pd<-input_pd[!is.na(input_pd$target_group),]
  message(print(summary(as.factor(input_pd$target_group))))
  ####################################

  # define the break time and colors
  if(break_month == "auto"){
    break_month<-break_month(input = input_pd$time,block = 6)
  }
  max_month<- break_month*6
  ###########################################
  if(is.null(cols))  cols<-  palettes(category = "box",palette = palette, show_col = F)


  # print(unique(input_pd$target_group))

  if(length(unique(input_pd$target_group))>2){
    input_pd<-input_pd[!is.na(input_pd$target_group),]
    # input_pd$target_group<-ifelse(input_pd$target_group=="High",3,ifelse(input_pd$target_group=="Middle",2,1))
    # Sur <- Surv(input_pd$time,input_pd$status)
    # pvalue<-getHRandCIfromCoxph(coxph(Surv(input_pd$time,input_pd$status)~input_pd$target_group,data = input_pd))
    # HR <- paste("Hazard Ratio = ", round(pvalue[,2],2), sep = "")
    # CI <- paste("95% CI: ", paste(round(pvalue[,3],2), round(pvalue[,4],2), sep = " - "), sep = "")
    ##########################################
    sfit <-  survminer:: surv_fit(Surv(input_pd$time,input_pd$status) ~ input_pd$target_group,data=input_pd)

    # hack strata for better survival curve
    names(sfit$strata) <- gsub("target_group=", "", names(sfit$strata))
    #############################################
    #############################################
    pp<- survminer:: ggsurvplot(sfit,
                                data             = input_pd,
                                censor           = TRUE,
                                ncensor.plot     = F,conf.int = F,
                                xlim             = c(0,max_month),
                                break.time.by    = break_month,
                                xlab             = "Months after diagnosis",
                                # legend.labs    = c(paste0('Low ',mini_sig),paste0("Middle ",mini_sig),paste0("High ",mini_sig)),
                                submain          = paste0(target_group," ",project),
                                surv.median.line = "h", # draw horizontal line at median survival
                                risk.table       = T,
                                tables.height    = 0.25,
                                palette          = cols,
                                pval.size        = 8)

    fitd <- survdiff(Surv(time,status)~ target_group,
                     data= input_pd,
                     na.action = na.exclude)

    p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

    # add nominal pvalue for log-rank test
    p.lab <- paste0("Overall P",
                    ifelse(p.val < 0.001, " < 0.001",
                           paste0(" = ",round(p.val, 3))))

    pp$plot <- pp$plot + annotate("text",
                                  x = 0, y = 0.55,
                                  hjust = 0,
                                  fontface = 3,
                                  # size = 1,
                                  label = p.lab)

    # calculate pair-wise survival comparison
    ps <- pairwise_survdiff(Surv(time,status)~ target_group,
                            data= input_pd,
                            p.adjust.method = "none")

    # add pair-wise comparison table
    # options(stringsAsFactors = FALSE)
    addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                             round(ps$p.value, 3))))
    addTab[is.na(addTab)] <- "-"
    # options(stringsAsFactors = TRUE)

    df <- tibble(x = 0, y = 0, tb = list(addTab))
    pp$plot <- pp$plot + ggpp::geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE, size = font.size.table)

  }else{

    # Sur <-   Surv(time = input_pd$time,event = input_pd$status)
    # turn high to '1'

    # print(unique(input_pd$target_group))
    levels <- unique(input_pd$target_group)
    print(levels)
    #####################################
    if(! reference_group%in% c(input_pd$target_group)) stop(">>>--- Reference_group must be one of target...")

    input_pd$target_group<-ifelse(input_pd$target_group == reference_group, 1, 0)

    pvalue<-getHRandCIfromCoxph(coxph(Surv(time = input_pd$time, event = input_pd$status)~input_pd$target_group,data = input_pd))

    HR <- paste("Hazard Ratio = ", round(pvalue[,2],2), sep = "")

    CI <- paste("95% CI: ", paste(round(pvalue[,3],2), round(pvalue[,4],2), sep = " - "), sep = "")
    # cut_off<-paste("cutoff = ",round(res.cut,3),sep = "")
    ###########################################
    sfit <- surv_fit(Surv(time = input_pd$time, event = input_pd$status) ~ input_pd$target_group,data = input_pd)
    # input_pd[,index]<-ifelse(input_pd[,index]==1,"High","LOW")
    ############################################
    # define break time


    levels<-levels[order(levels)]
    ###########################################
    pp<-survminer:: ggsurvplot(sfit,
                               data = input_pd,
                               censor = TRUE,
                               ncensor.plot = F,conf.int = F,
                               xlim = c(0,max_month),
                               break.time.by = break_month,
                               xlab = "Months after diagnosis",
                               legend.labs = c(paste0(levels[2]), paste0(levels[1])),
                               submain= paste0(target_group," ",project),
                               risk.table = T,
                               tables.height = 0.20,
                               palette = cols,
                               pval.size = 8,
                               pval = paste(pval = ifelse(pvalue[,1] < 0.0001, "P < 0.0001",
                                                          paste("P = ",round(pvalue[,1],4), sep = "")),
                                            HR, CI,sep = "\n"))

  }
  pp<-list(pp)
  res <- arrange_ggsurvplots(pp, print = FALSE, ncol = 1, nrow = 1)

  ggsave(res,filename = paste0(index,"-KMplot-",target_group,"-",project,".", fig.type),
         width = width ,height = height,path = save_path)
  return(res)
}





#' Break Time Period into Equal Blocks
#'
#' This function divides a given maximum time into specified equal blocks, allowing adjustments based on the time unit (e.g., months or days). It is useful for creating time intervals for further analysis or visualization.
#'
#' @param input Numeric vector or single numeric value representing time durations. 
#' @param block The number of blocks to divide the time into; default is 6 blocks.
#' @param time_type The unit of time in the input; either "month" for months or "day" for days, with default being "month".
#'
#' @return A numeric value representing the size of each block in the time unit specified.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Example with time in months
#' break_month(input = c(12, 24, 36), block = 6)
#'
#' # Example with time in days, converting to months
#' break_month(input = c(360, 720, 1080), block = 6, time_type = "day")
break_month<-function(input, block = 6, time_type = "month"){

  max_time <-max(input)
  if(time_type=="month"){
    max_time<-max_time
  }else if(time_type=="day"){
    max_time<-max_time/30
  }

  message(paste0("  Maximum of follow up time is ", max_time, " months; and will be divided into ", block, " sections;"))
  round(c(max_time %/% block)/5,0)*5
}






