





#' Title sig_box_batch
#' 
#' @description The sig_box_batch function is used to generate multiple signature box plots for analyzing the relationship between variables and groups in a given input data. It takes the input data, a vector of features or variables to analyze, and a vector of groups or categories to group the analysis by. The function allows customization of various parameters such as the path to save the figures, the appearance of the plots, and the statistical details to include. The generated figures are saved in the specified path or a default folder.
#'
#' @param input  The input data for analysis.
#' @param vars A vector of features or variables to analyze.
#' @param groups A vector of groups or categories to group the analysis by.
#' @param path The path to save the generated figures. If not specified, a default folder named "1-sig-box-batch" will be created.
#' @param index The index number to start naming the generated figures.
#' @param angle_x_text The angle of the x-axis labels in the figures.
#' @param hjust The horizontal alignment of the x-axis labels in the figures.
#' @param palette The color palette to use for the figures.
#' @param cols A vector of colors to use for the elements in the figures.
#' @param jitter A logical value indicating whether to add jitter to the points in the figures.
#' @param point_size The size of the points in the figures.
#' @param size_of_font The size of the font in the figures.
#' @param size_of_pvalue The size of the p-value text in the figures.
#' @param show_pvalue A logical value indicating whether to show the p-value in the figures.
#' @param return_stat_res A logical value indicating whether to return the statistical results.
#' @param assay The assay or data type to use for analysis (for Seurat object)
#' @param slot The slot in the data to use for analysis. (for Seurat object)
#' @param scale A logical value indicating whether to scale the data before analysis.
#' @param height The height of the generated figures.
#' @param width The width of the generated figures.
#' @param fig_type The file format to save the generated figures (e.g., "pdf").
#' @param pattern_vars A logical value indicating whether the feature names specified in 'vars' should be treated as patterns to match. If set to TRUE, the function will search for matching column names in the input data and use them for analysis.
#'
#' @return
#' @export
#'
#' @examples
sig_box_batch <- function(input, vars, groups, pattern_vars = FALSE, path = NULL, index = 0, angle_x_text = 0,
                          hjust = 0.5, palette = "jama", cols = NULL, jitter = FALSE, point_size = 5, size_of_font = 8,
                          size_of_pvalue = 4.5, show_pvalue = TRUE, return_stat_res = FALSE, assay = NULL, slot = "scale.data",
                          scale = FALSE, height = 5, width = 3.5, fig_type = "pdf", max_count_feas = 30){
  
  

  if(pattern_vars){
    vars_com <- c(NULL)
    for (i in 1:length(vars)) {
      vars_loop <- colnames(input)[str_detect(colnames(input), pattern = vars[i])]
      vars_com <- c(vars_com, vars_loop)
    }
    vars <- unique(vars_com[!is.na(vars_com)])
    if(length(vars)>max_count_feas) vars <- vars[1:max_count_feas]
    message(">>>== Variables that will be analysised :")
    message(paste0(vars, collapse = ", "))
  }
  ########################################
  if(length(vars)==1 & length(groups)==1) stop(">>>== `sig_box_batch` is suitable for cases where the `vars` or `groups` is greater than one ")
  
  if(is.null(path)){
    path <- creat_folder("1-sig-box-batch")
  }else{
    path <- creat_folder(path)
  }
  ########################################
  
  if(length(vars)>1){
    for (i in 1:length(vars)) {
      
      message(paste0(">>>== Processing feature: ", vars[i], "/n"))
      p<- sig_box(
          data            = input,
          signature       = vars[i],
          variable        = groups,
          angle_x_text    = angle_x_text,
          hjust           = hjust,
          palette         = palette,
          cols            = cols,
          jitter          = jitter,
          point_size      = point_size,
          size_of_font    = size_of_font,
          size_of_pvalue  = size_of_pvalue,
          show_pvalue     = show_pvalue,
          return_stat_res = return_stat_res,
          assay           = assay,
          slot            = slot,
          scale           = scale)
      
      ggsave(filename = paste0(index+i, "-", vars[i], "-and-", groups, ".", fig_type), height = height, width = width, path = path$folder_name)
    }
  }
  
  
  if(length(groups)>1){
    for (i in 1:length(groups)) {
      
      message(paste0(">>>== Processing feature: ", groups[i], "/n"))
      p<- sig_box(
        data            = input,
        signature       = vars,
        variable        = groups[i],
        angle_x_text    = angle_x_text,
        hjust           = hjust,
        palette         = palette,
        cols            = cols,
        jitter          = jitter,
        point_size      = point_size,
        size_of_font    = size_of_font,
        size_of_pvalue  = size_of_pvalue,
        show_pvalue     = show_pvalue,
        return_stat_res = return_stat_res,
        assay           = assay,
        slot            = slot,
        scale           = scale)
      
      height_index <- height
      
      levs <- unique(input[, groups[i]])
      levs <- length(levs[!is.na(levs)])
      height <- 2 + height_index * levs
      
      ggsave(filename = paste0(index+i, "-", vars, "-and-", groups[i], ".", fig_type), height = height, width = width, path = path$folder_name)
    }
  }
  
  message(paste0(">>>== Done"))
}
