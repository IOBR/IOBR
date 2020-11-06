

#' Transform signature gene data into list format
#'
#' @param sig_data signature data with signature name as column name and genes in row, using `NA` to substitude the missing value
#' @param save_signature logical, save siganture list as RData
#' @param output_name string of name to save signature RData
#'
#' @return list of signatures
#' @export
#' @author Dongqiang Zeng
#' @examples
format_signatures<-function(sig_data,save_signature = FALSE, output_name = "signatures"){

  message(paste0(">>> There are ",dim(sig_data)[2], "  signatures >>>" ))
  
  sig_data<-as.data.frame(sig_data)
  sig_data[sig_data=="NA"]<-NA
  #'reading each column and transfer it to a list
  bb<-as.list(NULL)
  for (i in 1:ncol(sig_data)) {
    aa<-as.character(sig_data[,i])
    aa<-list(aa);names(aa)<-names(sig_data[i])
    bb<-append(bb,aa)
  }
  bb<-lapply(bb,function(x) na.omit(x))
  bb<-lapply(bb,function(x) as.character(x))
  bb<-lapply(bb,function(x) unique(x))
  #' standerdized the name of list
  names(bb)<-gsub(names(bb),pattern = "\\.",replacement = "_")
  names(bb)<-gsub(names(bb),pattern = "\\ ",replacement = "_")
  names(bb)<-gsub(names(bb),pattern = "\\-",replacement = "_")
  # ########################################
  my_signatures<-bb
  ##########################################
  if(save_signature == TRUE){
    save(my_signatures,file = paste0(output_name,".RData"))
  }
  return(my_signatures)
}



