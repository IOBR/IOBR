




#' Calculate and Visualize Correlation between Two Variables
#' @description The "get_cor" function calculates and visualizes the correlation between two variables in a dataset. It provides options to scale the data, handle missing values, and incorporate additional data. The function supports various correlation methods and can display the correlation result. It generates a correlation plot with optional subtypes or categories, including a regression line. The plot can be customized with color palettes, labels, and titles. Additionally, the function allows saving the plot and data for further analysis.
#'
#' @param eset A dataset containing the variables of interest.
#' @param pdata An optional data frame providing additional data for the analysis. Default is NULL.
#' @param var1 The name of the first variable to analyze.
#' @param var2 The name of the second variable to analyze.
#' @param subtype An optional variable that represents subtypes/categories in the analysis. Default is NULL.
#' @param na.subtype.rm Whether to remove NA values in the subtype variable. Default is FALSE.
#' @param color_subtype The color palette for the subtypes/categories. Default is NULL.
#' @param palette The name of the color palette. Default is "jama".
#' @param index The index of the plot. Default is NULL.
#' @param method The correlation method to use. Default is "spearman".
#' @param show_cor_result Whether to print the correlation result. Default is TRUE.
#' @param col_line The color of the line in the correlation plot. Default is NULL.
#' @param id  The ID column name to use for labeling points on the plot. Default is "NULL".
#' @param show_lebel  Whether to show labels on the plot. Default is FALSE.
#' @param point_size The size of the points on the plot. Default is 4.
#' @param title The title of the plot. Default is NULL.
#' @param alpha The transparency level of points on the plot. Default is 0.7.
#' @param title_size The size of the title on the plot. Default is 2.
#' @param text_size The size of the text on the plot. Default is 15.
#' @param axis_angle The rotation angle of the axes labels on the plot. Default is 0.
#' @param hjust The horizontal justification of the axes labels on the plot. Default is 0.
#' @param show_plot Whether to display the plot. Default is TRUE.
#' @param fig.format The format of the saved figure. Default is "png".
#' @param fig.width The width of the saved figure. Default is 7.
#' @param fig.height The height of the saved figure. Default is 7.3.
#' @param path  The path where the figure will be saved. Default is NULL.
#' @param save_plot Whether to save the plot as a file. Default is FALSE.
#' @param id_pdata The column name in pdata that contains unique identifiers. Default is "ID".
#' @param scale Whether to scale the data. Default is TRUE.
#' @param is.matrix Whether the eset is a matrix data with feature as row names
#' @param id_eset The column name in eset that contains unique identifiers. Default is "ID".
#' @param add.hdr.line add.hdr.line
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data(eset_tme_stad, package = "IOBR")
#' get_cor(eset = eset_tme_stad, is.matrix = TRUE, var1 = "GZMB", var2 = "CD274")

get_cor <- function(eset, pdata = NULL, is.matrix = FALSE, id_eset = "ID", id_pdata = "ID", var1, var2, scale = TRUE,
        subtype = NULL, na.subtype.rm = FALSE, color_subtype = NULL,  palette = "jama", index = NULL,
        method = "spearman", show_cor_result = T, col_line = NULL, id = "NULL",
        show_lebel = FALSE, point_size = 4, title = NULL, alpha = 0.5,title_size = 1.5,
        text_size = 10, axis_angle = 0, hjust = 0, show_plot = TRUE, fig.format = "png",
        fig.width = 7, fig.height = 7.3, path = NULL, save_plot = FALSE, add.hdr.line = FALSE){

  if(is.null(index)) index<-1
  if(!show_lebel) id<-NULL
  #######################################
  if(is.null(pdata)){

    if(is.matrix){

      if(scale) eset <- scale_matrix(matrix = eset, log2matrix = FALSE, manipulate = FALSE)
      data<- as.data.frame(rownames_to_column(as.data.frame(t(eset)), var = "ID"))
      var1<- var1[var1%in%colnames(data)]
      var2<- var2[var2%in%colnames(data)]

    }else{
      colnames(eset)[which(colnames(eset)==id_eset)]<-"ID"
      data <- eset
      var1<- var1[var1%in%colnames(data)]
      var2<- var2[var2%in%colnames(data)]
      if(scale) data[,c(var1, var2)] <- scale(data[, c(var1, var2)])
    }

    # data[, unique(var1, var2)] <- apply(data[, unique(var1, var2)], 2, as.numeric)
    print(data[,c(var1, var2)])
  }else{
    colnames(pdata)[which(colnames(pdata)==id_pdata)]<-"ID"

    feas<- c(var1, var2)
    ##################################
    if(is.matrix){
      feas<- feas[feas%in%rownames(eset)]
      if(length(feas)==1) feas<- c(feas, "GAPDH", "TUBB1", "BMP2")

      data<- combine_pd_eset(eset = eset, pdata = pdata, id_pdata = "ID", feas = feas, scale = scale)

    }else{
      colnames(eset)[which(colnames(eset)==id_eset)]<-"ID"
      data <- eset
      data <- merge(pdata, eset, by ="ID")
    }

    message("Done: Combing data...")
    # print(head(data))
  }

  data<- data[!is.na(data[,var1]), ]
  data<- data[!is.na(data[,var2]), ]
  #######################################

  data<-as.data.frame(data)
  cor_result<- cor.test(data[,var1], data[,var2], method = method)
  if(show_cor_result) print(cor_result)

  pvalue<- exact_pvalue(data[,var1], data[,var2],method = method)
  print(paste0(">>>--- The exact p value is: ", pvalue))
  #######################################

  if(is.null(col_line)){

    if(cor_result$estimate>0){
      col_line<-"darkred"
    }else{
      col_line<-"steelblue"
    }
  }

  input2<-as.data.frame(data)

  if(is.null(color_subtype)){
    color_subtype<- palettes(category = "box",palette = palette, show_col = F, show_message = F)
  }

  if(!is.null(subtype)){
    colnames(input2)[which(colnames(input2)==subtype)]<-"categorys"

    if(na.subtype.rm) {
      input2<-input2[!is.na(input2$categorys),]
    }

    input2$categorys<- as.character(input2$categorys)
    input2[is.na(input2$categorys),"categorys"]<-"Not_avaliable"
    input2$categorys<-as.factor(input2$categorys)
    print(summary(input2$categorys))

    p<-ggplot(input2, aes(x = input2[,var1],
                          y =  input2[,var2],
                          colour = categorys)) +

      scale_color_manual(values= color_subtype)+
      geom_point(size = point_size,alpha = alpha) + # Show dots
      # scale_colour_gradient(high = "dark")+
      geom_smooth(method=lm, se=T,color= col_line)+
      labs(x = paste0(var1),
           y = paste0(var2),
           title = title,
           subtitle = paste0("r = ", round(unique(cor_result$estimate),3),
                             ",  P = ", sprintf("%1.1e", pvalue)))

  }


  if(is.null(subtype)){
    p<-ggplot(input2, aes(x = input2[,var1],
                          y =  input2[,var2])) +
      geom_point(size = point_size, alpha = alpha, colour = "black") + # Show dots
      # scale_colour_gradient(high = "dark")+
      geom_smooth(method=lm, se=T,color= col_line)+
      labs(x = paste0(var1),
           y = paste0(var2),
           title = title,
           subtitle = paste0("r = ", round(unique(cor_result$estimate),3),
                             ",  P = ", sprintf("%1.1e", pvalue)))
  }

  if(is.null(id) & show_lebel) stop("If you want to show label, parameter 'id' should be defined")

  if(!is.null(id)&show_lebel){
    p<-p+geom_text(label= as.character(data[, id]),
                   nudge_x = 0.25, nudge_y = 0.25, check_overlap = T)
  }

  theme<-design_mytheme(axis_title_size = title_size,
                        axis_text_size = text_size,
                        hjust = hjust,
                        axis_angle = axis_angle)

  p<-p+theme+theme(plot.subtitle=element_text(size=15, hjust= 0, face="italic", color="black"))

  if(add.hdr.line) {
    p <- p +  geom_hdr_lines()
  }

  if(show_plot) print(p)

  if(is.null(path)){
    path<-paste0("1-Cor-of-",var1,"-and-",var2)
  }

  ff <- creat_folder(path)

  if(save_plot){
    if(fig.format=="pdf"){
      ggsave(p,filename = paste0(index, "-",var2,"-",var1,"-correlation",".pdf"),
             width = fig.width ,height =fig.height,path = path)
    }else if(fig.format=="png"){
      ggsave(p,filename = paste0(index, "-",var2,"-",var1,"-correlation",".png"),
             width = fig.width ,height =fig.height,path = path)
    }

   save(data, file = paste0(ff$abspath, "0-input-data-",var1, "-", var2, ".RData"))
  }

  return(p)

}
