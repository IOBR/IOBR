





#' batch_surv - Batch Survival Analysis Function
#'
#' @description  This function is used to perform batch survival analysis. It calculates hazard ratios and confidence intervals for the specified variables based on the given data containing time-related information.
#'
#' @param pdata The data frame containing the data.
#' @param time  The column name representing time. Default is "time".
#' @param status The column name representing status. Default is "status".
#' @param variable The variables to perform survival analysis on.
#' @param best_cutoff Whether to use the best cutoff for survival analysis. Default is FALSE. If set to TRUE, the function will calculate hazard ratios and confidence intervals for the binary version of the variables using the best cutoff method.
#'
#' @author Dongqiang Zeng
#' @importFrom  survival coxph
#' @importFrom survival Surv
#' @export
#' @return A data frame containing hazard ratios and confidence intervals.
#' @examples
#'
#' # Loading TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Identifying signatures associated with gastric cancer patient survival
#' batch_surv(pdata = sig_stad, variable = colnames(sig_stad)[69:ncol(sig_stad)], time = "OS_time", status = "OS_status")

batch_surv <- function(pdata, variable, time ="time", status ="status", best_cutoff = FALSE){

  pdata<-as.data.frame(pdata)

  colnames(pdata)[which(colnames(pdata)==time)]<-"time"
  colnames(pdata)[which(colnames(pdata)==status)]<-"status"

  pdata<-pdata[!is.na(pdata$time),]
  pdata<-pdata[!is.na(pdata$status),]

  pdata$time <-as.numeric(pdata$time)
  pdata$status <-as.numeric(pdata$status)
  #################################################
  if(best_cutoff==TRUE){
   for (i in 1:length(variable)) {
     var <- variable[i]
     pdata <- best_cutoff(pdata = pdata, time = "time", status = "status", variable = var, PrintResult = FALSE)
     pdata[, paste0(var, "_binary")] <- ifelse(pdata[, paste0(var, "_binary")]=="High", 1, 0)
   }
    # print(pdata)
    variable <- paste0(variable, "_binary")
  }
  #################################################
  result_list<-lapply(pdata[,variable], function(x) coxph(Surv(pdata$time,pdata$status)~x,data = pdata[,variable]))
  result<-data.frame(NULL)
  for (i in 1:length(result_list)) {
    result1<-getHRandCIfromCoxph(result_list[[i]])
    rownames(result1)<-variable[i]
    result<-rbind(result,result1)
  }
  # result[result>1000]<-Inf

  result<-result[order(result$P,decreasing = F),]
  result<-tibble::rownames_to_column(result,var = "ID")

  if(best_cutoff==TRUE){
    result$ID <- gsub(result$ID, pattern = "_binary", replacement = "")
  }
  result<-tibble::as_tibble(result)
  return(result)
}


