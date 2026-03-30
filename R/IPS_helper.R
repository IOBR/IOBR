#' Calculate Immunophenoscore
#'
#' @description
#' Maps input score to Immunophenoscore (IPS) on a 0-10 scale.
#'
#' @param x A numeric value representing the input score.
#'
#' @return A numeric value between 0 and 10.
#'
#' @export
#'
#' @examples
#' ips <- ipsmap(2.5)
#' ips <- ipsmap(-1)
#' ips <- ipsmap(5)
ipsmap <- function(x) {
  if (is.na(x) || x <= 0) {
    0
  } else if (x >= 3) {
    10
  } else {
    round(x * 10 / 3, digits = 0)
  }
}


#' Map Score to Color
#'
#' @description
#' Maps a numeric input value to a color from a predefined palette.
#'
#' @param x A numeric value to be mapped to a color.
#' @param my_palette Color palette vector (should have 1001 colors).
#'   Default uses red-white-blue gradient.
#'
#' @return A color from the palette.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Create palette with 1001 colors
#' my_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(1001)
#' color <- mapcolors(2, my_palette)
#' color <- mapcolors(-2, my_palette)
#' }
mapcolors <- function(x, my_palette = NULL) {
  if (is.null(my_palette)) {
    my_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(1001)
  }

  if (x >= 3) {
    za <- 1000
  } else if (x <= -3) {
    za <- 1
  } else {
    za <- round(166.5 * x + 500.5, digits = 0)
  }

  my_palette[za]
}


#' Map Score to Black and White Color
#'
#' @description
#' Maps a numeric input value to a color from a black and white palette.
#'
#' @param x A numeric value to be mapped to a color.
#' @param my_palette2 Color palette vector (should have 1001 colors).
#'   Default uses black-white gradient.
#'
#' @return A color from the black and white palette.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Create palette with 1001 colors
#' my_palette2 <- grDevices::colorRampPalette(c("black", "white"))(1001)
#' color <- mapbw(1.5, my_palette2)
#' color <- mapbw(-1, my_palette2)
#' }
mapbw <- function(x, my_palette2 = NULL) {
  if (is.null(my_palette2)) {
    my_palette2 <- grDevices::colorRampPalette(c("black", "white"))(1001)
  }

  if (x >= 2) {
    za2 <- 1000
  } else if (x <= -2) {
    za2 <- 1
  } else {
    za2 <- round(249.75 * x + 500.5, digits = 0)
  }

  my_palette2[za2]
}
