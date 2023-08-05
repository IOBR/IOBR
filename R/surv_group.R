





#' Drawing KM-plot to a Categorical variable with two or more groups
#'
#' @param input_pdata DATA
#' @param project name of main title
#' @param time column name of follow up time
#' @param status status
#' @param time_type month or day
#' @param break_month break of time
#' @param palette default is `jama`, if cols is null
#' @param save_path default is KMplot
#' @param mini_sig prefix of variable
#' @param target_group target of binary variable
#' @param levels level of variable
#' @param reference_group default is
#' @param index folder name
#' @param cols if null, palette should be set
#' @param fig.type default is pdf, other option is png
#' @param ID identifier of data
#' @param width default is 6
#' @param height default is 6.5
#' @param font.size.table default is 62
#'
#' @author Dongqiang Zeng
#' @import survival
#' @return
#' @export
#'
#' @examples
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
  if(time!="time" & "time"%in%colnames(input_pdata)) stop(paste0("time already exists in column name of input_pdata."))
  if(status!="status" & "status"%in%colnames(input_pdata)) stop(paste0("status already exists in column name of input_pdata."))

  #filter data
  colnames(input_pdata)[which(colnames(input_pdata)==time)]<-"time"
  colnames(input_pdata)[which(colnames(input_pdata)==status)]<-"status"
  colnames(input_pdata)[which(colnames(input_pdata)==ID)]<-"ID"


  input_pdata$time<-as.numeric(input_pdata$time)
  input_pdata$status<-as.numeric(as.character(input_pdata$status))
  input_pdata<- input_pdata %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::filter(!is.na(status)) %>%
    dplyr::filter(!.$time=="NA") %>%
    dplyr::filter(!.$status=="NA") %>%
    dplyr::filter(.$time > 0)
  ###################################
  input_pdata$time<-as.numeric(input_pdata$time)
  input_pdata$status<-as.numeric(input_pdata$status)
  ###################################
  input_pdata<-as.data.frame(input_pdata[,colnames(input_pdata)%in%c("ID","time","status",target_group)])
  input_pdata$ID<-as.character(input_pdata$ID)
  # input_pdata$ProjectID<-project
  ##################################

  #transform follow up time
  if(time_type=="day") input_pdata$time<-input_pdata$time/30
  ###################################
  #cut use quantile
  #################################
  message(paste(">>> Dataset's survival follow up time is range between",
                paste(round(summary(input_pdata$time),2)[c(1,6)],collapse = " to "),"months"))


  ###################################
  save(input_pdata,file = paste0(abspath,"0-",project,"-",target_group,"-survival-analysis-input.RData"))
  #######################################

  colnames(input_pdata)[which(colnames(input_pdata)==target_group)]<-"target_group"
  input_pdata<-input_pdata[!is.na(input_pdata$target_group),]
  message(print(summary(as.factor(input_pdata$target_group))))
  ####################################

  # define the break time and colors
  if(break_month == "auto"){
    break_month<-break_month(input = input_pdata$time,block = 6)
  }
  max_month<- break_month*6
  ###########################################
  if(is.null(cols))  cols<-  palettes(category = "box",palette = palette, show_col = F)


  # print(unique(input_pdata$target_group))

  if(length(unique(input_pdata$target_group))>2){
    input_pdata<-input_pdata[!is.na(input_pdata$target_group),]
    # input_pdata$target_group<-ifelse(input_pdata$target_group=="High",3,ifelse(input_pdata$target_group=="Middle",2,1))
    # Sur <- Surv(input_pdata$time,input_pdata$status)
    # pvalue<-getHRandCIfromCoxph(coxph(Surv(input_pdata$time,input_pdata$status)~input_pdata$target_group,data = input_pdata))
    # HR <- paste("Hazard Ratio = ", round(pvalue[,2],2), sep = "")
    # CI <- paste("95% CI: ", paste(round(pvalue[,3],2), round(pvalue[,4],2), sep = " - "), sep = "")
    ##########################################
    sfit <-  survminer:: surv_fit(Surv(input_pdata$time,input_pdata$status) ~ input_pdata$target_group,data=input_pdata)

    # hack strata for better survival curve
    names(sfit$strata) <- gsub("target_group=", "", names(sfit$strata))
    #############################################
    #############################################
    pp<- survminer:: ggsurvplot(sfit,
                                data             = input_pdata,
                                censor           = TRUE,
                                ncensor.plot     = F,conf.int = F,
                                xlim             = c(0,max_month),
                                break.time.by    = break_month,
                                xlab             = "Months after diagnosis",
                                # legend.labs    = c(paste0('Low ',mini_sig),paste0("Middle ",mini_sig),paste0("High ",mini_sig)),
                                submain          = paste0(target_group,"-in-",project),
                                surv.median.line = "h", # draw horizontal line at median survival
                                risk.table       = T,
                                tables.height    = 0.25,
                                palette          = cols,
                                pval.size        = 8)

    fitd <- survdiff(Surv(time,status)~ target_group,
                     data= input_pdata,
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
                            data= input_pdata,
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

    # Sur <-   Surv(time = input_pdata$time,event = input_pdata$status)
    # turn high to '1'

    # print(unique(input_pdata$target_group))
    levels <- unique(input_pdata$target_group)
    print(levels)
    #####################################
    if(! reference_group%in% c(input_pdata$target_group)) stop(">>>--- Reference_group must be one of target...")

    input_pdata$target_group<-ifelse(input_pdata$target_group == reference_group, 1, 0)

    pvalue<-getHRandCIfromCoxph(coxph(Surv(time = input_pdata$time, event = input_pdata$status)~input_pdata$target_group,data = input_pdata))

    HR <- paste("Hazard Ratio = ", round(pvalue[,2],2), sep = "")

    CI <- paste("95% CI: ", paste(round(pvalue[,3],2), round(pvalue[,4],2), sep = " - "), sep = "")
    # cut_off<-paste("cutoff = ",round(res.cut,3),sep = "")
    ###########################################
    sfit <- surv_fit(Surv(time = input_pdata$time, event = input_pdata$status) ~ input_pdata$target_group,data = input_pdata)
    # input_pdata[,index]<-ifelse(input_pdata[,index]==1,"High","LOW")
    ############################################
    #' define break time


    levels<-levels[order(levels)]
    ###########################################
    pp<-survminer:: ggsurvplot(sfit,
                               data = input_pdata,
                               censor = TRUE,
                               ncensor.plot = F,conf.int = F,
                               xlim = c(0,max_month),
                               break.time.by = break_month,
                               xlab = "Months after diagnosis",
                               legend.labs = c(paste0(levels[2]), paste0(levels[1])),
                               submain= paste0(target_group,"-in-",project),
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





#' break month
#'
#' @param input
#' @param block
#' @param time_type
#'
#' @author Dongqiang Zeng
#' @return
#' @export
#'
#' @examples
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






