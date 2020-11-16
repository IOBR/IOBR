





##' Batch survival analysis on a data frame with multiple variables
##'
##'  @description  Extract hazard ratio and confidence intervals from a coxph object of subgroup analysis
##'
##' @param pdata data with follow up data and multiple variables
##' @param time  survival time
##' @param status survival status
##' @param variable use column c(a:b) to select variables
##'
##' @author Dongqiang Zeng
##' @importFrom  survival coxph
##' @importFrom survival Surv
##' @export
##' @examples
##'
##' sig_surv_result<- batch_surv(pdata = pdata_sig_tme_binary,variable <- c(c(ncol(pdata_stad)+1):ncol(pdata_sig_tme_binary)))

batch_surv<-function(pdata,variable,time="time", status="status"){

  pdata<-as.data.frame(pdata)
  colnames(pdata)[which(colnames(pdata)==time)]<-"time"
  colnames(pdata)[which(colnames(pdata)==status)]<-"status"

  result_list<-lapply(pdata[,variable], function(x) coxph(Surv(pdata$time,pdata$status)~x,data = pdata[,variable]))

  result<-data.frame(NULL)
  for (i in 1:length(result_list)) {
    result1<-getHRandCIfromCoxph(result_list[[i]])
    rownames(result1)<-colnames(pdata)[variable[i]]
    result<-rbind(result,result1)
  }
  result[result>100]<-Inf
  result<-result[order(result$P,decreasing = F),]
  result<-tibble::rownames_to_column(result,var = "ID")
  result<-tibble::as_tibble(result)
  return(result)
}


