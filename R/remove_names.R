


#' @export
patterns_to_na<-c("_cibersort","xCell","_EPIC","_TIMER","_quantiseq","_MCP","HALLMARK_", "_CIBERSORT","xcell","_timer","_mcp", "_epic")



#' Batch to transform patterns to specical character
#'
#' @param input_df input data frame
#' @param variables colnames or names of variable
#' @param patterns_to_na patterns that will be transform into NA: ''
#' @param patterns_space patterns that will be transform into space: ' '
#'
#' @return
#' @export
#'
#' @examples
remove_names<-function(input_df, variable = "colnames", patterns_to_na = patterns_to_na, patterns_space = NULL){


  if(variable=="colnames"){
    for (i in 1:length(patterns_to_na)) {
      colnames(input_df)<-gsub(colnames(input_df),pattern = patterns_to_na[i],replacement = "")
    }

    if(!is.null(patterns_space)){
      for (j in 1:length(patterns_space)) {
        colnames(input_df)<-gsub(colnames(input_df),pattern = patterns_space[j],replacement = " ")
      }
    }

  }else{
    input_df<-as.data.frame(input_df)
    for (i in 1:length(patterns_to_na)) {
      input_df[,variable]<-gsub(input_df[,variable],pattern = patterns_to_na[i],replacement = "")
    }
    if(!is.null(patterns_space)){
      for (j in 1:length(patterns_space)) {
        input_df[,variable]<-gsub(input_df[,variable],pattern = patterns_space[j],replacement = " ")
      }
    }
  }
  return(input_df)
}
