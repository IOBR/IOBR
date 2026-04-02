#' Map Score to Immunophenoscore
#'
#' @description
#' Maps input score to Immunophenoscore (IPS) on a 0-10 scale. Scores ≤0 map to 0,
#' scores ≥3 map to 10, and intermediate scores are linearly scaled.
#'
#' @param x Numeric value representing the aggregate z-score.
#'
#' @return Integer value between 0 and 10 representing the Immunophenoscore.
#'
#' @export
#'
#' @examples
#' ips <- ipsmap(2.5)
#' ips <- ipsmap(-1)
#' ips <- ipsmap(5)
ipsmap <- function(x) {
  if (is.na(x) || x <= 0) {
    0L
  } else if (x >= 3) {
    10L
  } else {
    as.integer(round(x * 10 / 3, digits = 0))
  }
}


#' Map Score to Color
#'
#' @description
#' Maps a numeric input value to a color from a blue-white-red gradient palette.
#' Values are mapped to a 1001-color palette where -3 maps to blue, 0 maps to white,
#' and +3 maps to red.
#'
#' @param x Numeric value to be mapped to a color (typically between -3 and 3).
#' @param my_palette Color palette vector (should have 1001 colors).
#'   Default uses blue-white-red gradient.
#'
#' @return A color from the palette as a hex code.
#'
#' @export
#'
#' @examples
#' my_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(1001)
#' color <- mapcolors(2, my_palette)
#' color <- mapcolors(-2, my_palette)
mapcolors <- function(x, my_palette = NULL) {
  if (is.null(my_palette)) {
    my_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(1001)
  }

  za <- if (x >= 3) {
    1000L
  } else if (x <= -3) {
    1L
  } else {
    as.integer(round(166.5 * x + 500.5, digits = 0))
  }

  my_palette[za]
}


#' Map Score to Black and White Color
#'
#' @description
#' Maps a numeric input value to a color from a black-white gradient palette.
#' Values are mapped to a 1001-color palette where -2 maps to black and +2 maps to white.
#'
#' @param x Numeric value to be mapped to a color (typically between -2 and 2).
#' @param my_palette2 Color palette vector (should have 1001 colors).
#'   Default uses black-white gradient.
#'
#' @return A color from the black-white palette as a hex code.
#'
#' @export
#'
#' @examples
#' my_palette2 <- grDevices::colorRampPalette(c("black", "white"))(1001)
#' color <- mapbw(1.5, my_palette2)
#' color <- mapbw(-1, my_palette2)
mapbw <- function(x, my_palette2 = NULL) {
  if (is.null(my_palette2)) {
    my_palette2 <- grDevices::colorRampPalette(c("black", "white"))(1001)
  }

  za2 <- if (x >= 2) {
    1000L
  } else if (x <= -2) {
    1L
  } else {
    as.integer(round(249.75 * x + 500.5, digits = 0))
  }

  my_palette2[za2]
}
