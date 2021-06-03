





#' Forest plot of survival analysis results.
#'
#' @param signature column of signatures 
#' @param data survival result
#' @param pvalue P value of survival result
#' @param HR Hazard ratio of survival result
#' @param CI_low_0.95 95% confidence interval 
#' @param CI_up_0.95 95% confidence interval 
#' @param n maximum of signatures
#' @param max_character maximum of character
#' @param discrete_width discrete width, default is 35
#' @param color_option default is 1, other option:2, 3
#'
#' @return
#' @export
#'
#' @examples
#' sig_surv_result<- batch_surv(pdata = pdata_sig_tme_binary,variable <- c(100:ncol(pdata_sig_tme_binary)))
#' sig_forest(data = sig_surv_result, signature = "ID")
sig_forest<-function(data, signature, pvalue = "P", HR = "HR", CI_low_0.95 = "CI_low_0.95", CI_up_0.95 = "CI_up_0.95", n = 10, max_character = 25, discrete_width = 35, color_option = 1) {
  
  
  data<-as.data.frame(data)
  colnames(data)[which(colnames(data)==signature)]<-"signature"
  colnames(data)[which(colnames(data)==pvalue)]<-"P"
  colnames(data)[which(colnames(data)==HR)]<-"HR"
  colnames(data)[which(colnames(data)==CI_low_0.95)]<-"CI_low_0.95"
  colnames(data)[which(colnames(data)==CI_up_0.95)]<-"CI_up_0.95"
  
  
  data[,c("P","HR","CI_low_0.95","CI_up_0.95")]<-apply(data[,c("P","HR","CI_low_0.95","CI_up_0.95")],2,function(x) as.numeric(x))
  
  data<- data[complete.cases(data),]
  #######################################
  
  if(dim(data)[1]>n){
    message(paste0("Top ", n," signatures will be shown"))
    
    data$HR_statistic<-data$HR - 1
    good_features<-high_var_fea(result = data,
                                target = "signature",
                                name_padj = "P",
                                padj_cutoff = 1,
                                name_logfc = "HR_statistic",
                                logfc_cutoff  = 0,
                                n = n/2)
    data<-data[data$signature%in%good_features,]
    data<-data[order(data$HR,decreasing = T),]
    data$signature<-as.character(data$signature)
  }
  #######################################
  
  if(max(nchar(data$signature))> max_character){
    
    goi<-as.character(data$signature)
    
    for (i in 1:length(goi)) {
      
      if(nchar(goi[i])>max_character){
        data[i, "signature"]<- gsub(data[i, "signature"],pattern = "\\_",replacement = " ")
      }
      
    }
  }
  
  
  
  pp<-ggplot(data=data,aes(x = HR,y = signature, color = P))+
    geom_errorbarh(aes(xmax= CI_up_0.95,xmin = CI_low_0.95),color="black",height=0,size=1.2)+
    geom_point(aes(x = HR,y= signature),size=6,shape=16)+
    geom_vline(xintercept = 1,linetype='dashed',size=0.8)+
    scale_x_continuous(breaks = c(0.5,1,1.50))+
    coord_trans(x='log2')+
    ylab("Signatures")+
    xlab(paste0("Hazard ratios of signatures"))+
    labs(color="P value")+
    viridis::scale_color_viridis(option= color_option)+
    
    theme_light()+
    theme(axis.text.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+
    theme(axis.text.y = element_text(size = 16, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+
    theme(title = element_text(size = 15,colour = "black", vjust = 0.5, hjust = 0.5))+
    scale_y_discrete(labels=function(x) stringr::str_wrap(x, width = discrete_width))
  
 print(pp)
 return(pp)
}
