





##' Batch survival analysis on a data frame with multiple variables
##'
##'  @description  Extract hazard ratio and confidence intervals from a coxph object of subgroup analysis
##'
##' @param pdata data with follow up data and multiple variables
##' @param time  survival time
##' @param status survival status
##' @param variable column names of included variables
##' @param best_cutoff
##'
##' @author Dongqiang Zeng
##' @importFrom  survival coxph
##' @importFrom survival Surv
##' @export
##' @examples
##'
##' sig_surv_result<- batch_surv(pdata = pdata_sig_tme_binary,variable <- c(c(ncol(pdata_stad)+1):ncol(pdata_sig_tme_binary)))

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
     # pdata[, var] <- ifelse(pdata[, var]=="High", 1, 0)
   }
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


