




#' save signature data into a data frame
#'
#' @param signatures output signatures, a list or a string
#' @param format output format
#' @param file.name file name of output data
#'
#' @return
#' @export
#'
#' @examples

output_sig<-function(signatures, format = "csv", file.name ){
  
  
  if(length(signatures)<=1){
    sig<-as.data.frame(sig)
    colnames(sig)<-"signature genes"
  }
  
  if(length(signatures)>=2){
    sig<-as.data.frame(t(do.call("rbind", signature_tme)))
    
    for (i in 1:ncol(sig)) {
      
      index<-duplicated(sig[,i])
      sig[index,i]<-NA
      
    }
  }
  
  if(format=="csv") write.csv(sig, file = paste0(file.name,".csv"))
  
  if(format=="rdata") save(sig, file = paste0(file.name,".RData"))
  return(sig)
}

