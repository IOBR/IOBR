




#' Reading signature gene matrix as list
#'
#' @param file_name the file name of csv file with signature name as column name and genes in row
#' @param save_signature logical, save siganture list as RData
#' @param output_name string of name to save signature RData
#'
#' @return list of signatures
#' @export
#' @author Dongqiang Zeng
#' @examples
#' mysignature<-read_signatures(file_name = "signatures.csv",save_signature = TRUE,output_name = "mysignatures")

read_signatures<-function(file_name,save_signature = FALSE, output_name = "signatures"){


  #creat-signature-file-before using this function: csv file using NA to substitde the missing value
  mysig<-read.csv(paste0(file_name),header = T)
  message(paste0(">>> There are ",dim(mysig)[2], "  signatures >>>" ))

  #'reading each column and transfer it to a list
  bb<-as.list(NULL);bb
  for (i in 1:ncol(mysig)) {
    aa<-as.character(mysig[,i])
    aa<-list(aa);names(aa)<-names(mysig[i])
    bb<-append(bb,aa)
  }
  bb<-lapply(bb,function(x) na.omit(x))
  bb<-lapply(bb,function(x) as.character(x))
  bb<-lapply(bb,function(x) unique(x))
  #' standerdized the name of list
  names(bb)<-gsub(names(bb),pattern = "\\.",replacement = "_")
  # ########################################
  my_signatures<-bb
  print(my_signatures[1:3])
  ##########################################
  if(save_signature == TRUE){
    save(my_signatures,file = paste0(output_name,".RData"))
  }
  return(my_signatures)
}



