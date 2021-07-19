




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
#' @param show_heatmap_col_name  show heatmpa column name
#' @param path save path, folder name
#' @param index default is null
#'
#' @return
#' @export
#'
#' @examples
sig_heatmap<-function(input, ID = "ID", features, group, palette = 2, show_col = F, 
                      show_palettes = F, show_plot = T, width = 8,
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
  
  colnames(input)[which(colnames(input)==ID)]<-"ID"
  # input<-column_to_rownames(input,var = "ID")
  
  input<-input[,c("ID", group, features)]
  colnames(input)[which(colnames(input)==group)]<-"target_group"
  
  pf_long_group <- tidyr::pivot_longer(input, 3:ncol(input), names_to = "variables",values_to = "value")
  
  ###################################################
  pf_long_group$value[pf_long_group$value > 2.5] = 2.5
  pf_long_group$value[pf_long_group$value < -2.5] = -2.5
  
  height_heatmap<-length(features)*0.2 + 3
  ####################################################
  heatmap_col<-palettes(category = "tidyheatmap",palette = palette, show_col = show_col, show_message = show_palettes)
  ####################################################
  pp<-pf_long_group %>%
    group_by(target_group) %>%
    tidyHeatmap:: heatmap(
      .column = ID,
      .row = variables,
      .value = value,
      # column_title = group_name,
      # annotation = group2,
      palette_value = heatmap_col,
      show_column_names = show_heatmap_col_name)
  
  if(show_plot) print(pp)
  
  if(is.null(index)) index<-1
  pp %>%  tidyHeatmap::save_pdf(paste0(abspath, index, "-",group,"-tidyheatmap.pdf"),
                                width = width,
                                height = height_heatmap)
}



