

#' Assimilate data with similar columns
#'
#' @param data_a originial dataframe
#' @param data_b manipulat dataframe
#'
#' @return  data_b
#' @export
#'
#' @examples
#'
#' pdata_a<-as.data.frame(array(NA,c(5,9)));pdata_a
#' colnames(pdata_a)<-LETTERS[1:dim(pdata_a)[2]];pdata_a
#' pdata_b<-as.data.frame(array(NA,c(7,4)));pdata_b
#' colnames(pdata_b)<-c("A","C","E","F");pdata_b
#' assimilate_data(data_a = pdata_a,data_b = pdata_b)
#'
assimilate_data<-function(data_a,data_b){
  data_a<-as.data.frame(data_a)
  data_b<-as.data.frame(data_b)
  miss<-as.character(colnames(data_a)[!colnames(data_a)%in%colnames(data_b)])
  ncolumns<-length(miss)
  nrows<-dim(data_b)[1]
  missdata<-as.data.frame(array(NA,c(nrows,length(miss))))
  colnames(missdata)<-miss
  data_b<-as.data.frame(cbind(data_b,missdata))
  data_b<-data_b[,match(colnames(data_a),colnames(data_b))]
  return(data_b)
}



