





#' Title sig_pheatmap
#' @description The sig_pheatmap function is used to generate heatmaps. It takes a data frame as input and creates a heatmap with grouping variables based on the provided parameters. The parameters include the input data, the features to include in the heatmap, the grouping variable, and optional additional grouping variables. It also offers flexible options to specify colors, adjust the size and layout of the heatmap, and save output files. The function returns a list containing annotation data, cluster colors, the plot object, and the transformed input matrix, allowing users to further analyze and visualize the heatmap data.
#'
#' @param input This parameter represents the input data for the heatmap. It should be a data frame with variables in column
#' @param feas This parameter specifies the features to include in the heatmap. It should be a vector containing the names of the columns in the "input" data frame.
#' @param group This parameter indicates the grouping variable for the heatmap. It should be a column name in the "input" data frame.
#' @param group2  (optional): This parameter represents an additional grouping variable for the heatmap. It should be a column name in the "input" data frame.
#' @param group3  (optional): This parameter indicates another additional grouping variable for the heatmap. It should be a column name in the "input" data frame.
#' @param ID (default: "ID"): This parameter allows you to specify the column name to be used as the sample identifier. The default value is "ID".
#' @param path (optional): This parameter represents the path where the output files will be stored. If not provided, the default path will be used.
#' @param cols1  (default: "random"): This parameter specifies the colors to be used for the first grouping variable. It can either be a vector of color names or the value "random" to randomly generate colors. The default value is "random".
#' @param seed (default: 54321): This parameter specifies the seed for the random number generator used in color selection. The default value is 54321.
#' @param show_col (default: FALSE): This parameter determines whether to display the colors used in the heatmap. The default value is FALSE.
#' @param cluster_cols  (default: TRUE): This parameter determines whether to cluster the columns of the heatmap. The default value is TRUE.
#' @param palette_for_heatmape (default: 6): This parameter represents the palette number for the heatmap. The default value is 6.
#' @param scale.matrix  (default: TRUE): This parameter specifies whether to scale the input matrix. The default value is TRUE.
#' @param cellwidth (default: 9): This parameter determines the width (in points) of each cell in the heatmap. The default value is 9.
#' @param cellheight  (default: 9): This parameter specifies the height (in points) of each cell in the heatmap. The default value is 9.
#' @param fig.type (default: "pdf"): This parameter indicates the file format for saving the heatmap. The default value is "pdf".
#' @param width (default: 9): This parameter represents the width (in inches) of the saved heatmap. The default value is 9.
#' @param height (optional): This parameter specifies the height (in inches) of the saved heatmap. If not provided, a default height will be calculated based on the number of features.
#' @param file_name_prefix  (default: 1): This parameter allows you to specify a prefix for the saved file name. The default value is 1.
#' @param cols2 (default: "random"): This parameter represents the colors to be used for the second grouping variable. It can either be a vector of color names or the value "random" to randomly generate colors. The default value is "random".
#' @param cols3 (default: "random"): This parameter indicates the colors to be used for the third grouping variable. It can either be a vector of color names or the value "random" to randomly generate colors. The default value is "random".
#' @param palette1 (default: 1): This parameter represents the palette number for grouping variable 1. The default value is 1.
#' @param palette2 (default: 2): This parameter indicates the palette number for grouping variable 2. The default value is 2.
#' @param palette3 (default: 3): This parameter specifies the palette number for grouping variable 3. The default value is 3.
#' @param show_colnames (default: FALSE): This parameter determines whether to display colum names
#'
#' @return
#' @export
#'
#' @author Dongqiang Zeng
#' @examples
#' data("tcga_stad_sig", package = "IOBR")
#' data("tcga_stad_pdata", package = "IOBR")
#' input <- merge(tcga_stad_pdata, tcga_stad_sig, by = "ID")
#' feas <- colnames(input)[grep(colnames(input), pattern = "MCPcounter")]
#' sig_pheatmap(input = input, feas = feas , group = "subtype", scale.matrix = TRUE)
#'
sig_pheatmap<- function(input, feas, group,
             group2               = NULL,
             group3               = NULL,
             ID                   =  "ID",
             path                 = NULL,
             cols1                = "random",
             cols2                = "random",
             cols3                = "random",
             seed                 = 54321,
             show_col             = FALSE,
             palette1              = 1,
             palette2              = 2,
             palette3              = 3,
             cluster_cols         = TRUE,
             palette_for_heatmape = 6,
             scale.matrix         = TRUE,
             cellwidth            = 1,
             cellheight           = 9,
             show_colnames        = FALSE,
             fig.type             = "pdf",
             width                = 6,
             height               = NULL,
             file_name_prefix     = 1){


  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("Marker-heatmap-average")
  }

  path<- creat_folder(file_store)
  ###################################

  input<-as.data.frame(input)
  feas<-feas[feas%in%colnames(input)]

  colnames(input)[which(colnames(input)==ID)]<-"idd"
  ###################################

  eset<-input[,c("idd", feas)]
  rownames(eset)<-NULL
  eset<-column_to_rownames(eset, var = "idd")
  if(scale.matrix == TRUE) eset<- scale(eset, scale = T, center = T)
  eset<-t(eset)
  ##################################

 if(is.null(group3)&is.null(group2)){
   anno<- input[,c("idd", group)]
 }else if(!is.null(group2)&is.null(goup3)){
   anno<- input[,c("idd", group, group2)]
 }else{
   anno<- input[,c("idd", group, group2, group3)]
 }

 rownames(anno) <- anno$idd
 anno$idd <- NULL
 ###################################


  mapal <- palettes(category = "heatmap", palette = palette_for_heatmape, counts = 200,show_col = show_col)

  if(is.null(height)){
    # if(is.null(group)) stop("group must be define")
    height<- 2 + length(feas)*0.25
  }
  ####################################################


  ##############################################
  if(length(cols1)==1){
    if(cols1=="random"){

      mycols1<-palettes(category = "random", palette = palette1, show_col = show_col, show_message = FALSE)
      message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
      set.seed(seed)
      mycols1<-mycols1[sample(length(mycols1), length(mycols1))]
      if(show_col) scales::show_col(mycols1)

    }else if(cols1 == "normal"){
      mycols1<-palettes(category = "random", palette = palette1, show_col = show_col, show_message = FALSE)
    }
  }else{
    mycols1<-cols1
    if(show_col) scales::show_col(mycols1)
  }
  ####################################################

  ##############################################
  if(!is.null(group2)){
    if(length(cols2)==1){
      if(cols2=="random"){

        mycols2<-palettes(category = "random", palette = palette2, show_col = show_col, show_message = FALSE)
        message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
        set.seed(seed)
        mycols2<-mycols2[sample(length(mycols2), length(mycols2))]
        if(show_col) scales::show_col(mycols2)

      }else if(cols2 == "normal"){
        mycols2<-palettes(category = "random", palette = palette2, show_col = show_col, show_message = FALSE)
      }
    }else{
      mycols2<-cols2
      if(show_col) scales::show_col(mycols2)
    }

  }

  ##############################################
  if(!is.null(group3)){
    if(length(cols3)==1){
      if(cols3=="random"){

        mycols3<-palettes(category = "random", palette = palette3, show_col = show_col, show_message = FALSE)
        message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
        set.seed(seed)
        mycols3<-mycols3[sample(length(mycols3), length(mycols3))]
        if(show_col) scales::show_col(mycols3)

      }else if(cols3 == "normal"){
        mycols3<-palettes(category = "random", palette = palette3, show_col = show_col, show_message = FALSE)
      }
    }else{
      mycols3<-cols3
      if(show_col) scales::show_col(mycols3)
    }
  }

  cluster_colors1<- mycols1[1:length(base:: unique(input[,group]))]
  names(cluster_colors1)<-unique(input[,group])
  cluster_colors1 <- list(group = cluster_colors1)
  names(cluster_colors1) <- group
  ######################################################
  if(!is.null(group2)){
    cluster_colors2<- mycols2[1:length(base:: unique(input[,group2]))]
    names(cluster_colors2)<-unique(input[,group2])
    cluster_colors2 <- list(group2 = cluster_colors2)
    names(cluster_colors2) <- group2
  }
  ######################################################
  if(!is.null(group3)){
    cluster_colors3<- mycols3[1:length(base:: unique(input[,group3]))]
    names(cluster_colors3)<- unique(input[,group3])
    cluster_colors3 <- list(group3 = cluster_colors3)
    names(cluster_colors3) <- group3
  }
  ######################################################

  if(is.null(group3)&is.null(group2)){
    cluster_colors<- cluster_colors1
  }else if(!is.null(group2)&is.null(goup3)){
    cluster_colors<- append(cluster_colors1, cluster_colors2)
  }else{
    cluster_colors<- append(cluster_colors1, cluster_colors2, cluster_colors3)
  }

  print(cluster_colors)
  #######################################################
  # library(pheatmap)
  # pdf(paste0(path$abspath, "2-markers-heatmap-of-average-", group, ".", fig.type), width = width, height = height)
  p<-pheatmap:: pheatmap(
    eset,
    color             = mapal, #colorRampPalette(c("darkblue", "white", "red3"))(99), #heatmap color
    scale             = "none",
    cluster_rows      = T,
    cluster_cols      = cluster_cols, #clustered by columns
    cellwidth         = cellwidth,
    cellheight        = cellheight,
    treeheight_col    = 6,
    treeheight_row    = 6,
    clustering_method = "complete",
    show_rownames     = T, #show cluster names
    show_colnames     = show_colnames,
    angle_col         = "45",
    # annotation_row    = annotation_row,
    annotation_col    = anno,
    annotation_colors = cluster_colors,
    fontsize          = 6,
    silent            =  T) #The color for clusters are sames as previous setting
  # dev.off()
  # print(p)
  p
  ggsave(p, filename =  paste0(file_name_prefix, "-pheatmap-", group, ".", fig.type),
         width = width,
         height = height,
         path = path$folder_name)

  res<-list( "p_anno" = anno, "p_cols" = cluster_colors,  "plot" = p, "eset" = eset)
  return(res)

}
