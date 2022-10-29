




#' Design my ggplpt2 theme
#'
#' @param plot_title_size title size
#' @param axis_title_size axis title size
#' @param axis_text_size axis text size
#' @param axis_angle  axis angle
#' @param legend.position legend.position default is bottom
#' @param legend.direction legend.direction horizontal or vertical
#' @param hjust horizontal justification
#' @param legend.size font size of legend
#' @param theme default is light
#' @param legend.size.text default is 12
#' @param legend.key.height default is 0.5
#' @param legend.key.width default is 0.5
#'
#' @author Dongqiang Zeng
#' @return mytheme
#' @export
#'
#' @examples
#' mytheme <- design_mytheme()
design_mytheme<-function(theme            = "light",
                         plot_title_size  = 2,
                         axis_title_size  = 2,
                         axis_text_size   = 12,
                         axis_angle       = 60,
                         hjust            = 1,
                         legend.position  = NULL,
                         legend.direction = NULL,
                         legend.size      = 0.25,
                         legend.key.height= 0.5,
                         legend.key.width = 0.5,
                         legend.size.text = 10){

  if(theme=="light"){
    mytheme<-ggplot2::theme_light()
  }else if(theme=="bw"){
    mytheme<-ggplot2:: theme_bw()
  }else if(theme=="classic"){
    mytheme<-ggplot2:: theme_classic()
  }else if(theme=="classic2"){
    mytheme<-ggpubr:: theme_classic2()
  }
  message(paste0(">>>>Options for `theme`: light, bw, classic and classic2"))

  mytheme<-mytheme + ggplot2:: theme(plot.title=element_text(size=rel(plot_title_size),hjust=0.5),
                      axis.title=element_text(size=rel(axis_title_size)),
                      axis.text=element_text(size=rel(4)),
                      axis.text.x= element_text(face="plain",size= axis_text_size,
                                                angle=axis_angle,hjust = hjust,color="black"),#family="Times New Roman"
                      axis.text.y= element_text(face="plain",size= axis_text_size,color="black"),#family="Times New Roman"
                      # panel.grid.major=element_line(color="white"),
                      # panel.grid.minor=element_line(color="white"),
                      # panel.border=element_rect(color="white"),
                      axis.line=element_line(color="black",size=0.5))

  if(is.null(legend.position)) {
    legend.position<- "bottom"
  }else{
    message(paste0(">>>>Options for 'legend.position' : none, left, right, bottom, top"))
  }

  if(is.null(legend.direction)) {
    legend.direction<- "horizontal"
  }else{
    message(paste0(">>>>Options for 'legend.direction' : horizontal, vertical "))
  }


  mytheme<-mytheme+
    theme(
      # legend.key.size=unit(legend.size,"cm"),
      legend.key.height= unit(legend.key.height, 'cm'),
      legend.key.width= unit(legend.key.width, 'cm'),
      # legend.title=element_blank(),
      legend.position= legend.position,#"none","left","right","bottom","top",or #c(0.5,1)
      legend.direction= legend.direction,# "vertical"
      legend.justification=c(.5,.5),#"center" or two-element numeric vector
      legend.box="horizontal",#"horizontal",
      legend.box.just="top",
      legend.text=element_text(colour="black", size= legend.size.text, face = "plain")#("plain", "italic", "bold", "bold.italic")
    )

  return(mytheme)
}



