#' Select Color Palettes for Visualization
#'
#' @description
#' Provides curated qualitative, sequential, and diverging palettes for multiple
#' plot types. Supports intensity adjustment and preview.
#'
#' @param category Character. Plot/palette category: one of `box`, `continue2`,
#'   `continue`, `random`, `heatmap`, `heatmap3`, `tidyheatmap`.
#' @param palette Character or numeric. Palette name or index (varies by category).
#' @param alpha Numeric. Alpha (transparency) scaling factor. Default is 1.
#' @param counts Integer. Number of colors (for continuous palettes).
#'   Default is 50.
#' @param show_col Logical. If TRUE, prints the palette. Default is TRUE.
#' @param show_message Logical. If TRUE, prints available options.
#'   Default is FALSE.
#'
#' @return Character vector of hex color codes.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' colors <- palettes(category = "box", palette = "nrc", show_col = FALSE)
#' heatmap_colors <- palettes(
#'   category = "heatmap", palette = 1, counts = 10, show_col = FALSE
#' )
palettes <- function(category = "box", palette = "nrc", alpha = 1,
                     counts = 50, show_col = TRUE, show_message = FALSE) {
  rlang::check_installed("ggsci")
  rlang::check_installed("RColorBrewer")
  rlang::check_installed("scales")

  valid_categories <- c(
    "box", "continue2", "continue", "random", "heatmap",
    "heatmap3", "tidyheatmap"
  )
  category <- rlang::arg_match(category, valid_categories)

  if (show_message) {
    cli::cli_alert_info(
      "Available categories: {paste(valid_categories, collapse = ', ')}"
    )
  }

  mypal <- switch(category,
    box = .palette_box(palette, alpha, show_message),
    continue2 = .palette_continue2(palette, alpha, show_message),
    random = .palette_random(palette),
    continue = .palette_continue(palette, show_message),
    heatmap = .palette_heatmap(palette, counts),
    heatmap3 = .palette_heatmap3(palette),
    tidyheatmap = .palette_tidyheatmap(palette, show_message)
  )

  if (show_col && !is.null(mypal)) {
    message(paste0("'", mypal, "'", collapse = ", "))
    scales::show_col(mypal)
  }

  mypal
}

#' @keywords internal
.palette_box <- function(palette, alpha, show_message) {
  valid_palettes <- c(
    "nrc", "jama", "aaas", "jco", "paired1", "paired2",
    "paired3", "paired4", "accent", "set2"
  )

  if (show_message) {
    cli::cli_alert_info(
      "Box palettes: {paste(valid_palettes, collapse = ', ')}"
    )
  }

  if (palette == "nrc") {
    return(ggsci::pal_npg("nrc", alpha = alpha)(9))
  }
  if (palette == "jama") {
    return(ggsci::pal_jama(palette = "default", alpha = alpha)(7))
  }
  if (palette == "aaas") {
    return(ggsci::pal_aaas(palette = "default", alpha = alpha)(9))
  }
  if (palette == "jco") {
    return(ggsci::pal_jco(palette = "default", alpha = alpha)(9))
  }
  if (palette == "paired1") {
    cols <- RColorBrewer::brewer.pal(11, "Paired")
    return(cols[1:8])
  }
  if (palette == "paired2") {
    cols <- RColorBrewer::brewer.pal(11, "Paired")
    return(cols[3:10])
  }
  if (palette == "paired3") {
    cols <- RColorBrewer::brewer.pal(11, "Paired")
    return(cols[5:11])
  }
  if (palette == "paired4") {
    cols <- RColorBrewer::brewer.pal(11, "Paired")
    return(cols[7:11])
  }
  if (palette == "accent") {
    return(RColorBrewer::brewer.pal(8, "Accent"))
  }
  if (palette == "set2") {
    return(RColorBrewer::brewer.pal(7, "Set2"))
  }

  cli::cli_abort("Palette {.val {palette}} not found. Use: {valid_palettes}")
}

#' @keywords internal
.palette_continue2 <- function(palette, alpha, show_message) {
  valid_palettes <- c("nrc", "jama", "aaas", "jco", "rdbu")

  if (show_message) {
    cli::cli_alert_info(
      "Continue2 palettes: {paste(valid_palettes, collapse = ', ')}"
    )
  }

  if (palette == "nrc") {
    cols <- ggsci::pal_npg("nrc", alpha = alpha)(9)
    return(cols[c(4, 1)])
  }
  if (palette == "jama") {
    cols <- ggsci::pal_jama(palette = "default", alpha = alpha)(9)
    return(cols[c(1, 4)])
  }
  if (palette == "aaas") {
    cols <- ggsci::pal_aaas(palette = "default", alpha = alpha)(9)
    return(cols[c(1, 6)])
  }
  if (palette == "jco") {
    cols <- ggsci::pal_jco(palette = "default", alpha = alpha)(9)
    return(cols[c(1, 2)])
  }
  if (palette == "rdbu") {
    cols <- RColorBrewer::brewer.pal(11, "RdBu")
    return(cols[c(10, 2)])
  }

  cli::cli_abort("Palette {.val {palette}} not found. Use: {valid_palettes}")
}

#' @keywords internal
.palette_random <- function(palette) {
  cli::cli_alert_info(
    "Random palettes: 1 (palette1), 2 (palette2), 3 (palette3), 4 (palette4)"
  )
  palette_num <- ifelse(palette %in% 1:4, palette, 4)
  load_data(paste0("palette", palette_num))
}

#' @keywords internal
.palette_continue <- function(palette, show_message) {
  valid_palettes <- c("rdbu", "puor", "blues", "reds")

  if (show_message) {
    cli::cli_alert_info(
      "Continue palettes: {paste(valid_palettes, collapse = ', ')}"
    )
  }

  palette_upper <- toupper(palette)

  if (palette == "rdbu") {
    return(RColorBrewer::brewer.pal(11, "RdBu"))
  }
  if (palette == "puor") {
    return(RColorBrewer::brewer.pal(11, "PuOr"))
  }
  if (palette == "blues") {
    return(RColorBrewer::brewer.pal(11, "Blues"))
  }
  if (palette == "reds") {
    return(RColorBrewer::brewer.pal(11, "Reds"))
  }

  cli::cli_abort("Palette {.val {palette}} not found. Use: {valid_palettes}")
}

#' @keywords internal
.palette_heatmap <- function(palette, counts) {
  cli::cli_alert_info(paste0(
    "Heatmap palettes: 1 (pheatmap), 2 (peach), 3 (blues), ",
    "4 (virids), 5 (reds), 6 (RdBu), 7 (navy_firebrick), 8 (magma)"
  ))

  if (palette == 1) {
    return(rev(colorRampPalette(
      RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")
    )(counts)))
  }
  if (palette == 2) {
    return(colorRampPalette(c("#3182bd", "white", "#dd1c77"))(counts))
  }
  if (palette == 3) {
    return(rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(counts)))
  }
  if (palette == 4) {
    return(grDevices::colorRampPalette(c(
      "#000004FF", "#51127CFF", "#B63679FF",
      "#FCA50AFF", "#F7F419FF"
    ))(counts))
  }
  if (palette == 5) {
    return(rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Reds"))(counts)))
  }
  if (palette == 6) {
    return(rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(counts)))
  }
  if (palette == 7) {
    return(colorRampPalette(c("navy", "white", "firebrick"))(counts))
  }
  if (palette == 8) {
    return(grDevices::colorRampPalette(c(
      "#000004FF", "#6B0A7DFF", "#B63679FF",
      "#FCA50AFF", "#F7F419FF"
    ))(counts))
  }

  cli::cli_abort("Heatmap palette {.val {palette}} must be 1-8")
}

#' @keywords internal
.palette_heatmap3 <- function(palette) {
  cli::cli_alert_info(paste0(
    "Heatmap3 palettes: pheatmap, peach, blues, virids, reds, normal"
  ))

  if (palette == "pheatmap") {
    return(c("#4575B4", "#FEF9B6", "#D73027"))
  }
  if (palette == "peach") {
    return(c("#3182bd", "white", "#dd1c77"))
  }
  if (palette == "blues") {
    return(c("#F7FBFF", "#88BEDC", "#084594"))
  }
  if (palette == "virids") {
    return(c("#000004FF", "#CD4347FF", "#FCFFA4FF"))
  }
  if (palette == "reds") {
    return(c("#FFF5F0", "#FB7555", "#99000D"))
  }
  if (palette == "normal") {
    return(c("navy", "white", "firebrick"))
  }

  cli::cli_abort(paste0(
    "Heatmap3 palette {.val {palette}} not found. ",
    "Use: pheatmap, peach, blues, virids, reds, normal"
  ))
}

#' @keywords internal
.palette_tidyheatmap <- function(palette, show_message) {
  if (show_message) {
    cli::cli_alert_info("Tidyheatmap palettes: 1-6")
  }

  if (palette == 1) {
    return(c("#000004FF", "#B63679FF", "#F7F419FF"))
  }
  if (palette == 2) {
    return(rev(RColorBrewer::brewer.pal(3, "RdYlBu")))
  }
  if (palette == 3) {
    return(rev(RColorBrewer::brewer.pal(3, "RdYlGn")))
  }
  if (palette == 4) {
    return(c("#9E0142", "#FFFFBF", "#5C4A42"))
  }
  if (palette == 5) {
    return(rev(RColorBrewer::brewer.pal(3, "PiYG")))
  }
  if (palette == 6) {
    return(c("navy", "white", "firebrick"))
  }

  cli::cli_abort("Tidyheatmap palette {.val {palette}} must be 1-6")
}
