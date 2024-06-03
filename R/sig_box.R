




#' Title: sig_box - Generate a Box plot with Statistical Comparisons
#'
#' @description The sig_box function is designed to generate a boxplot with optional statistical comparisons. It takes in various parameters such as data, signature, variable, and more to customize the plot. It can be used to visualize and analyze data in a Seurat object or any other data frame.
#' @param data The input data, which can be a Seurat object or a data frame.
#' @param signature The column name representing the y-axis values in the plot.
#' @param variable The column name representing the x-axis values in the plot.
#' @param palette The color palette to use for filling the boxplots (default = "nrc").
#' @param cols Optional vector of colors to use for filling the boxplots (default = NULL).
#' @param jitter Boolean flag indicating whether to add jitter to the data points (default = FALSE).
#' @param point_size The size of the data points when jitter is enabled (default = 5).
#' @param angle_x_text The angle (in degrees) at which x-axis text labels should be displayed (default = 0).
#' @param hjust The horizontal justification of x-axis text labels (default = 0.5).
#' @param show_pvalue Boolean flag indicating whether to display the statistical comparison results (default = TRUE).
#' @param return_stat_res Boolean flag indicating whether to return the statistical comparison results (default = FALSE).
#' @param size_of_pvalue The font size for the statistical comparison results (default = 6).
#' @param size_of_font The font size for axis labels and text (default = 10).
#' @param assay The name of the assay to extract data from in a Seurat object (default = NULL, uses the default assay).
#' @param slot The slot name to extract data from in a Seurat object (default = "scale.data").
#' @param scale Boolean flag indicating whether to scale the signature values (default = FALSE).
#'
#' @return Depending on `return_stat_res`, returns either a ggplot object of the boxplot or the results of the statistical comparison.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' data("tcga_stad_pdata", package = "IOBR")
#'
#' sig_box(data = tcga_stad_pdata, signature = "TMEscore_plus", variable = "subtype", jitter = T, palette = "jco")
#'
#' sig_box(data = tcga_stad_pdata, signature = "TMEscore_plus", variable = "subtype", jitter = FALSE, palette = "jco")
sig_box<-function(data, signature, variable, angle_x_text = 0, hjust = 0.5, palette = "nrc", cols = NULL, jitter = FALSE, point_size = 5, size_of_font = 10,
                  size_of_pvalue = 6, show_pvalue = TRUE, return_stat_res = FALSE, assay = NULL, slot = "scale.data", scale = FALSE){

  if(class(data)[1]=="Seurat"){

    cat(crayon::green(">>>-- Derive matrix data from Seurat object...\n"))
    input<- extract_sc_data( sce                = data,
                             vars               = signature,
                             assay              = assay,
                             slot               = slot,
                             combine_meta_data  = TRUE)
    data<- input
  }


  data<-as.data.frame(data)
  data<-data[,c(variable, signature)]
  colnames(data)[which(colnames(data)==variable)]<-"variable"
  colnames(data)[which(colnames(data)==signature)]<-"signature"

  if(scale) data[,"signature"] <-as.numeric(scale(data[,"signature"], scale = T, center = T))
  data<-data[!is.na(data$variable),]
  ####################################
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

  size_font <- size_of_font*0.2

  if(angle_x_text%in%c(30, 45, 60)) hjust = 1

  p<-p+theme_light()+
    theme(axis.title.y=element_text(size=rel(size_font)),
          axis.title.x = element_text(size=rel(size_font)),
          axis.text=element_text(size=rel(size_font)),
          axis.text.x= element_text(face="plain", angle = angle_x_text, hjust = hjust, color="black"), #family="Times New Roman"

          axis.line=element_line(color="grey", size=0.05))+
    theme(legend.position = "none")


  if(show_pvalue)   p<-p+stat_compare_means(comparisons = comparision,size = size_of_pvalue)+stat_compare_means(size = size_of_pvalue)


  res<-compare_means(signature ~ variable, data = data)
  print(res)

  if(jitter){
    p<- p + geom_point(shape = 21,
                       size = point_size,
                       position = position_jitterdodge(dodge.width = 0.2),
                       alpha = .5)
  }

  print(p)

  if(return_stat_res){
    return(res)
  }else{
    return(p)
  }

}
