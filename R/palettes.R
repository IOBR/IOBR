#' Select Color Palettes for Visualization
#'
#' Provides curated qualitative, sequential, and diverging palettes for multiple plot types. Supports intensity adjustment and preview.
#'
#' @param category Character. Plot/palette category: one of `box`, `continue2`, `continue`, `random`, `heatmap`, `heatmap3`, `tidyheatmap`.
#' @param palette Character or numeric. Palette name or index (varies by category).
#' @param alpha Numeric. Alpha (transparency) scaling factor. Default is 1.
#' @param counts Integer. Number of colors (for continuous palettes). Default is 50.
#' @param show_col Logical. If TRUE, prints the palette. Default is TRUE.
#' @param show_message Logical. If TRUE, prints available options. Default is FALSE.
#'
#' @return Character vector of hex color codes.
#' @author Dongqiang Zeng
#' @export
#' @examples
#' colors <- palettes(category = "box", palette = "nrc", show_message = TRUE)
#' heatmap_colors <- palettes(category = "heatmap", palette = 1, counts = 100, show_col = TRUE)
palettes <- function(category = "box", palette = "nrc", alpha = 1, counts = 50, show_col = TRUE, show_message = FALSE) {
  if (show_message) message(paste0("There are seven categories you can choose: box, continue2, continue, random, heatmap, heatmap3, tidyheatmap "))

  if (category == "box") {
    if (show_message) message(paste0("There are ten palettes you can choose: nrc, jama, aaas, jco, paired1-4, accent, set2"))
    if (palette == "nrc") {
      mypal <- ggsci::pal_npg("nrc", alpha = alpha)(9)
    } else if (palette == "jama") {
      mypal <- ggsci::pal_jama(palette = c("default"), alpha = alpha)(7)
    } else if (palette == "aaas") {
      mypal <- ggsci::pal_aaas(palette = c("default"), alpha = alpha)(9)
    } else if (palette == "jco") {
      mypal <- ggsci::pal_jco(palette = c("default"), alpha = alpha)(9)
    } else if (palette == "paired1") {
      mypal <- RColorBrewer::brewer.pal(11, "Paired")
      mypal <- mypal[1:8]
    } else if (palette == "paired2") {
      mypal <- RColorBrewer::brewer.pal(11, "Paired")
      mypal <- mypal[3:10]
    } else if (palette == "paired3") {
      mypal <- RColorBrewer::brewer.pal(11, "Paired")
      mypal <- mypal[5:11]
    } else if (palette == "paired4") {
      mypal <- RColorBrewer::brewer.pal(11, "Paired")
      mypal <- mypal[7:11]
    } else if (palette == "accent") {
      mypal <- RColorBrewer::brewer.pal(8, "Accent")
    } else if (palette == "set2") {
      mypal <- RColorBrewer::brewer.pal(7, "Set2")
    }

    if (show_col) {
      print(paste0("'", mypal, "'", collapse = ", "))
      scales::show_col(mypal)
    }
  }

  if (category == "continue2") {
    if (show_message) message(paste0("There are five palettes you can choose: nrc, jama, aaas, jco, rdbu"))
    if (palette == "nrc") {
      mypal <- ggsci::pal_npg("nrc", alpha = alpha)(9)
      mypal <- mypal[c(4, 1)]
    } else if (palette == "jama") {
      mypal <- ggsci::pal_jama(palette = c("default"), alpha = alpha)(9)
      mypal <- mypal[c(1, 4)]
    } else if (palette == "aaas") {
      mypal <- ggsci::pal_aaas(palette = c("default"), alpha = alpha)(9)
      mypal <- mypal[c(1, 6)]
    } else if (palette == "jco") {
      mypal <- ggsci::pal_jco(palette = c("default"), alpha = alpha)(9)
      mypal <- mypal[c(1, 2)]
    } else if (palette == "rdbu") {
      mypal <- RColorBrewer::brewer.pal(11, "RdBu")
      mypal <- mypal[c(10, 2)]
    }
    if (show_col) {
      print(paste0("'", mypal, "'", collapse = ", "))
      scales::show_col(mypal)
    }
  }

  if (category == "random") {
    message(">>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4")

    ###########################################################
    if (palette == 1) {
      data("palette1")
      mypal <- palette1
    } else if (palette == 2) {
      data("palette2")
      mypal <- palette2
    } else if (palette == 3) {
      data("palette3")
      mypal <- palette3
    } else if (palette == 4) {
      data("palette4")
      mypal <- palette4
    } else {
      data("palette4")
      mypal <- palette4
    }

    if (show_col) {
      print(paste0("'", mypal, "'", collapse = ", "))
      scales::show_col(mypal)
    }
  }

  if (category == "continue") {
    if (show_message) message(paste0("There are four palettes you can choose: rdbu, puor, blues, reds"))
    if (palette == "rdbu") {
      mypal <- RColorBrewer::brewer.pal(11, "RdBu")
    } else if (palette == "puor") {
      mypal <- RColorBrewer::brewer.pal(11, "PuOr")
    } else if (palette == "blues") {
      mypal <- RColorBrewer::brewer.pal(11, "Blues")
    } else if (palette == "reds") {
      mypal <- RColorBrewer::brewer.pal(11, "Reds")
    }
    if (show_col) {
      print(paste0("'", mypal, "'", collapse = ", "))
      scales::show_col(mypal)
    }
  }

  if (category == "heatmap") {
    message(paste0("There are five palettes you can choose: 1 = pheatmap,  2 = peach,  3 = blues, 4 = virids, 5 = reds, 6 = RdBu, 7 = navy_firebrick"))

    if (palette == 1) {
      mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))(counts))
    } else if (palette == 2) {
      mypal <- colorRampPalette(c("#3182bd", "white", "#dd1c77"))(counts)
    } else if (palette == 3) {
      mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(counts))
    } else if (palette == 4) {
      mypal <- viridis::inferno(counts)
    } else if (palette == 5) {
      mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Reds"))(counts))
    } else if (palette == 6) {
      mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(counts))
    } else if (palette == 7) {
      mypal <- colorRampPalette(c("navy", "white", "firebrick"))(counts)
    } else if (palette == 8) {
      mypal <- colorRampPalette(viridis::magma(5))(counts)
    }


    if (show_col) {
      print(paste0("'", mypal, "'", collapse = ", "))
      scales::show_col(mypal)
    }
  }

  if (category == "heatmap3") {
    message(paste0("There are six palettes you can choose: pheatmap, virids, blues, reds, peach, normal"))

    if (palette == "pheatmap") {
      mypal <- c("#4575B4", "#FEF9B6", "#D73027")
    } else if (palette == "peach") {
      mypal <- c("#3182bd", "white", "#dd1c77")
    } else if (palette == "blues") {
      mypal <- c("#F7FBFF", "88BEDC", "#084594")
    } else if (palette == "virids") {
      mypal <- c("#000004FF", "#CD4347FF", "#FCFFA4FF")
    } else if (palette == "reds") {
      mypal <- c("#FFF5F0", "FB7555", "#99000D")
    } else if (palette == "normal") {
      mypal <- c("navy", "white", "firebrick")
    }
    if (show_col) {
      print(paste0("'", mypal, "'", collapse = ", "))
      scales::show_col(mypal)
    }
  }

  if (category == "tidyheatmap") {
    if (show_message) message(paste0("There are six palettes you can choose: 1, 2, 3, 4, 5, 6"))
    if (palette == 1) {
      mypal <- circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), viridis::magma(5))
    } else if (palette == 2) {
      mypal <- circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), rev(RColorBrewer::brewer.pal(n = 5, name = "RdYlBu")))
    } else if (palette == 3) {
      mypal <- circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), rev(RColorBrewer::brewer.pal(n = 5, name = "RdYlGn")))
    } else if (palette == 4) {
      mypal <- circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), rev(RColorBrewer::brewer.pal(n = 5, name = "Spectral")))
    } else if (palette == 5) {
      mypal <- circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), rev(RColorBrewer::brewer.pal(n = 5, name = "PiYG")))
    } else if (palette == 6) {
      mypal <- circlize::colorRamp2(c(-3, 0, 3), c("navy", "white", "firebrick"))
    }
  }

  return(mypal)
}
