




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
    signatures<-as.data.frame(signatures)
    colnames(signatures)<-"signature genes"
  }

  if(length(signatures)>=2){
    signatures<-as.data.frame(t(do.call("rbind", signatures)))

    for (i in 1:ncol(signatures)) {

      index<-duplicated(signatures[,i])
      signatures[index,i]<-NA

    }
  }

  if(format=="csv") write.csv(signatures, file = paste0(file.name,".csv"))

  if(format=="rdata") save(signatures, file = paste0(file.name,".RData"))
  return(signatures)
}

