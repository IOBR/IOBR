




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
#'
#' @author Dongqiang Zeng
#' @return mytheme
#' @export
#'
#' @examples
#' mytheme <- design_mytheme()
design_mytheme<-function( plot_title_size = 2,axis_title_size = 2,
                          axis_text_size = 12, axis_angle = 60,
                          hjust = 1,legend.position = "bottom",
                          legend.direction="horizontal", legend.size = 12){

  mytheme<-ggplot2:: theme_light()+   ###theme_bw()
   ggplot2:: theme(plot.title=element_text(size=rel(plot_title_size),hjust=0.5),
          axis.title=element_text(size=rel(axis_title_size)),
          axis.text=element_text(size=rel(2.5)),
          axis.text.x= element_text(face="plain",size= axis_text_size,
                                    angle=axis_angle,hjust = hjust,color="black"),#family="Times New Roman"
          axis.text.y= element_text(face="plain",size= axis_text_size,color="black"),#family="Times New Roman"
          # panel.grid.major=element_line(color="white"),
          # panel.grid.minor=element_line(color="white"),
          # panel.border=element_rect(color="white"),
          axis.line=element_line(color="black",size=0.5))+theme(
            legend.key.size=unit(.3,"inches"),
            legend.title=element_blank(),
            legend.position= legend.position,#"none","left","right","bottom","top",or #c(0.5,1)
            legend.direction= legend.direction,# "vertical"
            legend.justification=c(.5,.5),#"center" or two-element numeric vector
            legend.box="horizontal",#"horizontal",
            legend.box.just="top",
            legend.text=element_text(colour="black",size=legend.size,face = "plain")#("plain", "italic", "bold", "bold.italic")
          )
  return(mytheme)
}



