
#' Batch Signature Survival Plot
#'
#' This function generates survival plots for multiple projects based on input data and signature scores.
#'
#' @param input_pdata This parameter represents the input data for the survival analysis. It is expected to be a data frame.
#' @param signature This parameter specifies the column name in input_pdata that represents the target signature for survival analysis.
#' @param ID This optional parameter represents the column name in input_pdata that contains the unique identifier for each data point. The default value is "ID".
#' @param column_of_Project A character string specifying the name of the column representing the project IDs. Default is "ProjectID".
#' @param project This optional parameter represents the project name. The default value is "KM".
#' @param time This optional parameter represents the column name in input_pdata that contains the time variable for survival analysis. The default value is "time".
#' @param status This optional parameter represents the column name in input_pdata that contains the survival status variable. The default value is "status".
#' @param time_type This optional parameter specifies the time unit used in the analysis. Default is "day".Options are "day" and "month".
#' @param break_month This optional parameter specifies the interval at which the time axis breaks should be made in the Kaplan-Meier plot. The default value is "auto".
#' @param palette This optional parameter allows the user to specify the palette type if cols is not provided. The default value is "jama".
#' @param save_path This optional parameter specifies the directory path for saving the plot. The default value is "KM-plot".
#' @param mini_sig This optional parameter represents the label for the mini signature in the legend. The default value is "score".
#' @param show_col logical variable, if TRUE, color will be print and show in the R studio
#' @param cols This optional parameter allows the user to specify the color palette for the plot. If not specified, a default palette will be used.
#' @param fig.type This optional parameter specifies the file format for saving the plot. The default value is "png".
#' @param project A character string specifying the project name. Default is NULL.
#' @param time A character string specifying the name of the column representing the time variable. Default is "time".
#' @param status A character string specifying the name of the column representing the status variable. Default is "status".
#' @param time_type A character string specifying the type of time variable. Default is "day".Options are "day" and "month".
#' @param break_month A character string specifying the break interval for the time variable. Default is "auto".
#' @param palette A character string specifying the color palette to be used for the plots. Default is "jama".
#' @param cols A vector specifying custom colors for the plots. Default is NULL.
#' @param mini_sig A character string specifying the minimum signature score. Default is "score".
#' @param save_path A character string specifying the path where the plots will be saved. Default is "Multiple-KM-plot".
#' @param show_col A logical value indicating whether to show colors in the plots. Default is TRUE.
#' @param fig.type A character string specifying the file type for saving the plots. Default is "pdf".
#'
#' @importFrom survminer surv_cutpoint
#' @import survival
#' @author Dongqiang Zeng
#' 
#' @return A data frame containing the combined data from all projects after survival analysis.
#' @export
#'
#' @examples
#' # Loading TCGA-STAD microenvironment signature data
#' data("sig_stad", package = "IOBR")
#' # Generating survival plots for multiple projects
#' sig_stad <- as.data.frame(sig_stad)
#' result <- batch_sig_surv_plot(input_pdata = sig_stad, signature = "T.cells.CD8", ID = "ID", column_of_Project = "ProjectID", project = NULL,
#'                               time = "OS_time", status = "OS_status", time_type = "month", break_month = "auto", palette = "jama",
#'                               cols = NULL, mini_sig = "score", save_path = "Multiple-KM-plot", show_col = TRUE, fig.type = "pdf")
#' print(result)
batch_sig_surv_plot<-function(input_pdata,signature, ID = "ID", column_of_Project="ProjectID", project = NULL, time = "time", status = "status", time_type = "day", break_month = "auto",
                                palette="jama", cols = NULL, mini_sig = "score", save_path = "Multiple-KM-plot", show_col = TRUE, fig.type = "pdf"){
  
  input_pdata <- as.data.frame(input_pdata)
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
