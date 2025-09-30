#' calculate Immunophenoscore
#' Title
#'
#' @param x A numeric value representing the input score.
#'
#' @return A numeric value representing the Immunophenoscore (IPS), scaled between 0 and 10.
#' @export
#'
#' @examples
#' score <- 2.5
#' ips <- ipsmap(score)
#' print(ips)
ipsmap <- function(x) {
  # if(is.na(x)) x<-0
  if (x <= 0) {
    ips <- 0
  } else {
    if (x >= 3) {
      ips <- 10
    } else {
      ips <- round(x * 10 / 3, digits = 0)
    }
  }
  return(ips)
}



#' Map Colors Based on Input Value
#'
#' This function maps a numeric input value to a color from a predefined palette.
#'
#' @param x A numeric value to be mapped to a color.
#'
#' @return A color from the predefined palette corresponding to the input value.
#' @export
#'
#' @examples
#' # Define a palette with 1001 colors
#' my_palette <- colorRampPalette(c("blue", "white", "red"))(1001)
#' # Map a value to a color
#' color <- mapcolors(2)
#' print(color)
mapcolors <- function(x) {
  za <- NULL
  if (x >= 3) {
    za <- 1000
  } else {
    if (x <= -3) {
      za <- 1
    } else {
      za <- round(166.5 * x + 500.5, digits = 0)
    }
  }
  return(my_palette[za])
}

#' Map Black and White Colors Based on Input Value
#'
#' This function maps a numeric input value to a color from a predefined black and white palette.
#'
#' @param x A numeric value to be mapped to a color.
#'
#' @return A color from the predefined black and white palette corresponding to the input value.
#' @export
#'
#' @examples
#' # Define a black and white palette with 1001 colors
#' my_palette2 <- colorRampPalette(c("black", "white"))(1001)
#' # Map a value to a color
#' color <- mapbw(1.5)
#' print(color)
mapbw <- function(x) {
  za2 <- NULL
  if (x >= 2) {
    za2 <- 1000
  } else {
    if (x <= -2) {
      za2 <- 1
    } else {
      za2 <- round(249.75 * x + 500.5, digits = 0)
    }
  }
  return(my_palette2[za2])
}
