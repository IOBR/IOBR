






#' signature box plot
#'
#' @param data input data
#' @param signature name of signature
#' @param variable name of category variable
#' @param palette default is nature palette
#' @param cols default is null
#' @param jitter default is FALSE
#' @param point_size size of point
#'
#' @return
#' @export
#'
#' @examples
#' sig_box(data = tcga_stad_pdata, signature = "TMEscore_plus", variable = "subtype",jitter = T)
sig_box<-function(data, signature, variable, palette = "nrc", cols = NULL, jitter = FALSE, point_size = 5){
  
  
  data<-as.data.frame(data)
  data<-data[,c(variable, signature)]
  
  colnames(data)[which(colnames(data)==variable)]<-"variable"
  
  colnames(data)[which(colnames(data)==signature)]<-"signature"
  
  data<-data[!is.na(data$variable),]
  
  if(is.null(cols)){
    cols<-IOBR::palettes(category = "box", palette = palette, show_message = FALSE, show_col = FALSE)
  }else{
    cols<-cols
  }
  
  p<-ggplot(data,aes(x = variable,y = signature,fill = variable))+
    geom_boxplot(notch = F,outlier.shape = NULL,outlier.size = 0)+
    scale_fill_manual(values= cols)+
    ylab(signature)+
    xlab(variable)
  
  comparision<-combn(unique(as.character(data$variable)), 2, simplify=F)
  
  p<-p+theme_light()+
    theme(axis.title.y=element_text(size=rel(2.0)),
          axis.title.x = element_text(size=rel(2.0)),
          axis.text=element_text(size=rel(2.0)),
          axis.line=element_line(color="black", size=0.5))+
    theme(legend.position = "none")
    
  p<-p+stat_compare_means(comparisons = comparision,size=6)+stat_compare_means(size=6)
  
  if(jitter){
    p<- p + geom_point(shape = 21, 
                       size = point_size, 
                       position = position_jitterdodge(dodge.width = 0.2), 
                       alpha = .5) 
  }
  
  print(p)
  return(p)
  
}

