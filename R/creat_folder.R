




#' Title
#'
#' @param f1
#' @param f2
#' @param return
#' @param f3
#'
#' @return
#' @export
#'
#' @examples
creat_folder<-function(f1, f2 = NULL, f3 = NULL, return = NULL){


  if(!is.null(f3)){
    path<-file.path(getwd(), f1, f2, f3)
    if(!dir.exists(path)) dir.create(file.path(getwd(), f1, f2, f3), recursive = TRUE)
  }else if(!is.null(f2)){
    path<-file.path(getwd(), f1, f2)
    if(!dir.exists(path)) dir.create(file.path(getwd(), f1, f2), recursive = TRUE)
  }else{
    path<-file.path(getwd(), f1)
    if(!dir.exists(path)) dir.create(file.path(getwd(), f1), recursive = TRUE)
  }


  if(is.null(return)){

    if(!is.null(f3)){
      res<-list("folder_name" = paste0(f1, "/", f2, "/", f3),
                "abspath" = paste0(file.path(getwd(), f1, f2, f3), "/"))
    }else if(!is.null(f2)){
      res<-list("folder_name" = paste0(f1, "/", f2),
                "abspath" = paste0(file.path(getwd(), f1, f2), "/"))
    }else{
      res<-list("folder_name" = f1,
                "abspath" = paste0(file.path(getwd(), f1), "/"))
    }

  }else{

    if(return == 1){
      res<-list("folder_name" = f1,
                "abspath" = paste0(file.path(getwd(), f1), "/"))
    }else if(return == 2){
      res<-list("folder_name" = paste0(f1, "/", f2),
                "abspath" = paste0(file.path(getwd(), f1, f2), "/"))
    }else if(return == 3){
      res<-list("folder_name" = paste0(f1, "/", f2, "/", f3),
                "abspath" = paste0(file.path(getwd(), f1, f2, f3), "/"))
    }

  }

  return(res)
}
