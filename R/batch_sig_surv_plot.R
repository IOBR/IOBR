




#' Title
#'
#' @param input_pdata 
#' @param signature 
#' @param ID 
#' @param column_of_Project 
#' @param time 
#' @param status 
#' @param time_type 
#' @param break_month 
#' @param palette 
#' @param mini_sig 
#' @param save_path 
#' @param show_col 
#' @param fig.type 
#' @param project 
#' @param cols 
#'
#' @return
#' @export
#'
#' @examples
batch_sig_surv_plot<-function(input_pdata,signature, ID = "ID", column_of_Project="ProjectID", project = NULL, time = "time", status = "status", time_type = "day", break_month = "auto",
                                palette="jama", cols = NULL, mini_sig = "score", save_path = "Multiple-KM-plot", show_col = TRUE, fig.type = "pdf"){
  
  
  goi<-as.character(levels(as.factor(as.character(input_pdata[,column_of_Project]))))
  
  # message(summary(as.factor(as.character(input_pdata[,column_of_Project]))))
  
  colnames(input_pdata)[which(colnames(input_pdata)==column_of_Project)]<-"ProjectID"
  
  input_pdata[,"ProjectID"]<-as.character(input_pdata[,"ProjectID"])
  ####################################
  
  input_pdata_com<-data.frame(NULL)
  
  for (i in 1:length(goi)) {
    
    var<-goi[i]
    message(">>> Preprocessing dataset: ",var)
    pd<-input_pdata[input_pdata$ProjectID%in%var,]
    
    pd<-sig_surv_plot(
      input_pdata = pd,
      signature   = signature,
      project     = var,
      ID          = ID,
      time        = time,
      status      = status,
      time_type   = time_type,
      break_month = break_month,
      cols        = cols,
      palette     = palette,
      show_col    =show_col,
      mini_sig    = mini_sig,
      fig.type    = fig.type,
      save_path   = save_path,
      index       = i)
    
    input_pdata_com<-rbind(input_pdata_com, pd$data)
    
  }
  
  return(input_pdata_com)
  
}
