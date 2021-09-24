





#' Extract best cut-off and add new binary object to data frame
#'
#' @param pdata phenotype data with surcival information and features
#' @param variable selected featrues that will be classified into two groups
#' @param time column name of survival time
#' @param status column name of event of follow up
#' @param PrintResult logit, print results of survival analysis before and after classifying.
#'
#' @return pdata with binary variables
#' @export
#'
#' @examples
#'
best_cutoff2<-function(pdata,variable,time = "time",status = "status",PrintResult = T){
  
  pdata<-as.data.frame(pdata)
  colnames(pdata)[which(colnames(pdata)==time)]<-"time"
  colnames(pdata)[which(colnames(pdata)==status)]<-"status"
  
  pdata<-pdata[!is.na(pdata$time),]
  pdata<-pdata[!is.na(pdata$status),]
  
  pdata$time<-as.numeric(pdata$time)
  pdata$status<-as.numeric(pdata$status)
  
  y<-Surv(pdata$time,pdata$status)
  iscutoff<-surv_cutpoint(pdata, time = "time",event = "status",variables =variable)
  aa<-paste("best cutoff = ",iscutoff$cutpoint$cutpoint)
  plot(iscutoff,variable,palette="npg")
  
  bb<-base:: summary(coxph(y ~ pdata[, which(colnames(pdata)==variable)],data = pdata))
  
  variable2<-paste(variable,"_binary",sep = "")
  pdata[,variable2]<-ifelse(pdata[,variable]< iscutoff$cutpoint$cutpoint,"Low","High")
  pdata[,variable2]<-as.factor(pdata[,variable2])
  cc<-summary(pdata[,variable2])
  dd<-base:: summary(coxph(y~pdata[,which(colnames(pdata)==variable2)],data =pdata))
  if(PrintResult) {
    print(list(best_cutoff = aa ,cox_continuous_object=bb,summary_binary_variable = cc,cox_binary_object=dd))
  }
  return(list(pdata = pdata, best_cutoff = iscutoff$cutpoint$cutpoint,cox_continuous_object=bb,summary_binary_variable = cc,cox_binary_object=dd))
  
}
