



#'  Visualize the results as a stacked bar chart with tidyverse/ggplot2.
#'
#' @param input input  result of `CIBERSORT` or `quantiseq`
#' @param id patient identifier of input
#' @param title plot title
#' @param legend.position legend position
#' @param coord_filp logical variables, `coord_filp` of ggplot2
#'
#' @return
#' @export
#'
#' @examples
#' res<-cell_bar_plot(input = cibersort[1:20,],id = "ID")
#'
cell_bar_plot<- function(input, id = "ID", title = "Cell Fraction", legend.position = "bottom", coord_filp = TRUE){

  input<-as.data.frame(input)
  colnames(input)[which(colnames(input)==id)]<-"ID"

  if("ProjectID"%in%colnames(input))  input<-input[,-which(colnames(input)=="ProjectID")]
  if("index"%in%colnames(input))  input<-input[,-which(colnames(input)=="index")]
  if("RMSE_CIBERSORT"%in%colnames(input))  input<-input[,-which(colnames(input)=="RMSE_CIBERSORT")]
  if("P-value_CIBERSORT"%in%colnames(input))  input<-input[,-which(colnames(input)=="P-value_CIBERSORT")]
  if("Correlation_CIBERSORT"%in%colnames(input)) input<-input[,-which(colnames(input)=="Correlation_CIBERSORT")]

  input<-remove_names(input_df = input, variable = "colnames", patterns_to_na = patterns_to_na, patterns_space = "_")
  ##################
  if(legend.position == "top"|legend.position=="bottom") {
    legend.direction<-"horizontal"
  }else{
    legend.direction<-"vertical"
  }

  if(coord_filp){
    pp<-input %>%
      gather(cell_type,fraction, -ID) %>%
      # plot as stacked bar chart
      ggplot(aes(x=ID, y=fraction, fill=cell_type)) +
      geom_bar(stat='identity') +
      coord_flip() +
      theme_light()+
      scale_fill_manual(values = palettes(category = "random",counts = 23, show_col = FALSE)) +
      scale_x_discrete(limits = rev(levels(input)))+
      ggtitle(paste0(title))+
      theme(plot.title=element_text(size=rel(2),hjust=0.5),
            axis.text.x= element_text(face="plain",angle=0,hjust = 1,color="black"),
            axis.text.y= element_text(face="plain",angle= 30,hjust = 1,color="black"))+
      theme(legend.title = element_blank(),
            legend.position= legend.position,
            legend.direction= legend.direction,
            legend.justification=c(.5,.5),
            legend.box="horizontal",
            legend.box.just="top")
  }else{
   pp<- input %>%
      gather(cell_type,fraction, -ID) %>%
      # plot as stacked bar chart
      ggplot(aes(x=ID, y=fraction, fill=cell_type)) +
      geom_bar(stat='identity') +
      # coord_flip() +
      theme_light()+
      scale_fill_manual(values = palettes(category = "random",counts = 23, show_col = FALSE)) +
      scale_x_discrete(limits = rev(levels(input)))+
      ggtitle(paste0(title))+
      theme(plot.title=element_text(size=rel(2),hjust=0.5),
            axis.text.x= element_text(face="plain",angle=0,hjust = 1,color="black"),
            axis.text.y= element_text(face="plain",angle=30,hjust = 1,color="black"))+
      theme(legend.title = element_blank(),
            legend.position= legend.position,
            legend.direction= legend.direction,
            legend.justification=c(.5,.5),
            legend.box="horizontal",
            legend.box.just="top")
  }

 print(pp)
 return(pp)

}



