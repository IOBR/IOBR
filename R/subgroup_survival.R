





##' Title Extract hazard ratio and confidence intervals from a coxph object of subgroup analysis
##'
##' @param pdata     data combine variables, follow up time, and outcome;
##' @param variables  object that have several levles;
##' @param time_name     column name of follow up time
##' @param status_name   column name of outcome event
##' @param object    target that use to calculate coxph fit;
##' If target has two levels which name `High` and `Low` group, please transform `High` into `0` and `Low` into `1` before executing the command
##' @author Dongqiang Zeng
##' @export
##' @import survival
##' @examples
##'
##' #source data and filter NA
##' data(subgroup_data)
##' input <- subgroup_data %>% filter(time > 0) %>% filter(!is.na(status)) %>% filter(!is.na(AJCC_stage))
##'
##' ##for binary variable
##' data1 <- subgroup_survival(pdata = input,time_name ="time", status_name = "status",
##' variables = c("ProjectID", "AJCC_stage"), object ="score_binary" )
##' data1
##'
##' ##for continue variables
##' data2 <- subgroup_survival(pdata = input, time_name = "time", status_name = "status",
##' variables = c("ProjectID", "AJCC_stage"), object = "score" )
##' data2

subgroup_survival<-function(pdata,time_name="time",status_name= "status",variables,object){
  P<-1;HR<-1;CI_low_0.95<-1;CI_up_0.95<-1;
  result<-data.frame(P,HR,CI_low_0.95,CI_up_0.95,row.names = "defult");
  pdata<-as.data.frame(pdata);
  for (sig in variables) {
    ind<-which(colnames(pdata)==sig);
    tmp<-pdata[!is.na(pdata[,ind]),];
    if(dim(tmp)[1]==0)
      next;
    tmp[,ind]<-as.factor(as.character(tmp[,ind]));

    result2<-data.frame(P,HR,CI_low_0.95,CI_up_0.95,row.names = "defult");
    nl<-nlevels(tmp[,ind]);
    for (i in 1:nl) {
      tmp1<-tmp[as.character(tmp[,ind])==names(summary(tmp[,ind]))[i],];
      if(dim(tmp1)[1]==0)
        next;
      fit<-coxph(Surv(tmp1[,time_name],tmp1[,status_name])~tmp1[,object],data = tmp1);
      result1<-getHRandCIfromCoxph(fit);
      rownames(result1)<-paste(sig,levels(tmp[,ind])[i],sep = "_");
      result2<-rbind(result2,result1);
    }
    result<-rbind(result,result2);
  }
  result[result>100]<-Inf
  return(result[-c(grep(rownames(result),pattern ="defult")),])
}
########################################
