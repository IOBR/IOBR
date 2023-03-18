




#' simple heatmap
#'
#' @param input data frame with `ID`, `group` and features
#' @param ID default is ID
#' @param features features
#' @param group category
#' @param palette default is 2, options: 1,2,3
#' @param show_col show color
#' @param show_palettes show palettes
#' @param show_plot show plot
#' @param width default is 8
#' @param show_heatmap_col_name  show heat map column name
#' @param path save path, folder name
#' @param index default is null
#' @param palette_group default is jama
#' @param height height of figure
#' @param size_col font size of column name
#' @param size_row font size of row name
#' @param angle_col angle of column name
#' @param cols_group user can define the colors of group
#' @param column_title  title of column
#' @param row_title title of row
#' @param scale defalut is FALSE
#'
#' @return
#' @export
#'
#' @examples
sig_heatmap<-function(input, ID = "ID", features, group, scale = FALSE, palette = 2, palette_group = "jama", show_col = F,
                      show_palettes = F, cols_group = NULL, show_plot = T, width = 8, height = NULL, size_col = 10,size_row = 8,
                      angle_col = 90, column_title = NULL, row_title = NULL,
                      show_heatmap_col_name = F, path = NULL, index = NULL){


  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("1-",group,"-relevant-varbiles-heatmap")
  }

  if(!file.exists(file_store)) dir.create(file_store)
  abspath<-paste(getwd(),"/",file_store,"/",sep ="" )


  input<-as.data.frame(input)
  features<-features[features%in%colnames(input)]

  colnames(input)[which(colnames(input)==ID)]<-"idd"
  # input<-column_to_rownames(input,var = "ID")

  input<-input[,c("idd", group, features)]
  colnames(input)[which(colnames(input)==group)]<-"target_group"

  pf_long_group <- tidyr::pivot_longer(input, 3:ncol(input), names_to = "variables",values_to = "value")

  ###################################################
  pf_long_group$value[pf_long_group$value > 2.5] = 2.5
  pf_long_group$value[pf_long_group$value < -2.5] = -2.5


  if(is.null(height)){
    height_heatmap<-length(features)*0.1 + 3
  }else{
    height_heatmap<- height
  }

  ####################################################
  heatmap_col<-palettes(category = "tidyheatmap", palette = palette, show_col = show_col, show_message = show_palettes)


  if(!is.null(cols_group)){
    color_box<-cols_group
  }else{
    if(is.null(palette)){
      color_box<-palettes(category = "random", show_col = show_col ,show_message = show_palettes)
    }else{
      color_box<-palettes(category = "box", palette = palette_group, show_col = show_col,show_message = show_palettes)
    }
  }

  target_level<- unique(as.character(pf_long_group$target_group))
  target_level<-target_level[!is.na(target_level)]

  n<- length(target_level)
  # print(n)
  color_box<-color_box[1:n]
  # print(color_box)
  ####################################################
 if(scale){
   scale = "column"
 }else{
   scale = "none"
 }

  pp<-pf_long_group %>%
    group_by(target_group) %>%
    tidyHeatmap:: heatmap(
      .column           = idd,
      .row              = variables,
      .value            = value,
      palette_grouping  = list(c(color_box)),
      scale             = scale,
      column_title      = column_title,
      row_title         = row_title,

      # annotation      = group2,
      palette_value     = heatmap_col,
      show_column_names = show_heatmap_col_name,
      column_names_gp   = grid::gpar(fontsize = size_col),
      row_names_gp      = grid::gpar(fontsize = size_row),
      column_names_rot  = angle_col)

  # if(!is.null(add_anno)){
  #   pp<-pp|>add_tile(add_anno, show_legend = TRUE)
  # }

  if(show_plot) print(pp)

  if(is.null(index)) index<-1
  pp %>%  tidyHeatmap::save_pdf(paste0(abspath, index, "-",group,"-tidyheatmap.pdf"),
                                width = width,
                                height = height_heatmap)

  return(pp)
}



