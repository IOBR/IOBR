




#' Merge Data Frames Handling Duplicated Column Names
#'
#' This function merges two data frames by resolving duplicated column names based on user preference. 
#' It provides options to choose which data frame's duplicated columns to keep, ensuring smooth merges 
#' without losing important data integrity.
#'
#' @param x The first data frame to be merged.
#' @param y The second data frame to be merged.
#' @param by.x The column name(s) in `x` used for merging.
#' @param by.y The column name(s) in `y` used for merging.
#' @param all.x Logical indicating if all rows from `x` should be included in the output.
#' @param all.y Logical indicating if all rows from `y` should be included in the output.
#' @param all Logical indicating if all rows from both `x` and `y` should be included in the output,
#'        superseding `all.x` and `all.y` if not NULL.
#' @param choose Specifies which data frame's duplicated non-joining columns should be retained: "x" or "y".
#'
#' @return A data frame resulting from merging `x` and `y` according to the specified parameters.
#' @export
#' @examples
#' df1 <- data.frame(ID = 1:3, Name = c("A", "B", "C"), Value = 1:3)
#' df2 <- data.frame(ID = 1:3, Name = c("X", "Y", "Z"), Score = 4:6)
#' merged_df <- merge_duplicate(df1, df2, by.x = "ID", by.y = "ID", all.x = TRUE, all.y = FALSE, choose = "x")
merge_duplicate<-function(x, y, by.x, by.y, all.x, all.y,all = NULL, choose = "x"){


  duplicate_names<- intersect(colnames(x),colnames(y))

  if(by.x == by.y)  duplicate_names<-duplicate_names[!duplicate_names%in%c(by.x)]

  if(choose == "x"){
    y<-y[,!colnames(y)%in%duplicate_names]
  }

  if(choose == "y"){
    x<-x[,!colnames(x)%in%duplicate_names]
  }
  ###################################
  if(!is.null(all)){
    res<-merge(x = x, y = y, by.x = by.x, by.y = by.y , all = all)
  }else{
    res<-merge(x = x, y = y, by.x = by.x, by.y = by.y , all.x = all.x, all.y = all.y)
  }
  return(res)

}
