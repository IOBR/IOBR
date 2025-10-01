#' Set and view the color palettes
#'
#' @param cols Vector of colors, users can define the cols manually.  This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param show_col Whether to show color palettes
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#'
#' @return Vector of colors
#' @export
#'
#' @examples
#' mycols <- get_cols()
get_cols <- function(cols = "normal", palette = 1, show_col = T, seed = 123) {
  ##############################################
  if (length(cols) == 1) {
    if (nchar(palette) > 1) {
      cols <- palettes(category = "box", palette = palette, show_col = show_col, show_message = FALSE)
      mycols <- cols
    } else {
      if (cols == "random") {
        mycols <- palettes(category = "random", palette = palette, show_col = show_col)
        message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
        set.seed(seed)
        mycols <- mycols[sample(length(mycols), length(mycols))]
        if (show_col) scales::show_col(mycols)
      } else if (cols == "normal") {
        mycols <- palettes(category = "random", palette = palette, show_col = show_col)
      }
    }
  } else {
    mycols <- cols
    if (show_col) scales::show_col(mycols)
  }

  return(mycols)
}
