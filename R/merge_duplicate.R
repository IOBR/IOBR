




#' Merge data frames with duplicated column names
#'
#' @param x
#' @param y
#' @param by.x
#' @param by.y
#' @param all.x
#' @param all.y
#' @param all
#' @param choose
#'
#' @return
#' @export
#'
#' @examples
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
