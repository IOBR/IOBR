#' Set and View Color Palettes
#'
#' @description
#' Retrieves color palettes from the IOBR package with options for randomization
#' and visualization. Users can specify predefined palettes or provide custom
#' colors.
#'
#' @param cols Character vector of colors, or one of:
#'   - `"normal"`: Use standard palette
#'   - `"random"`: Randomly shuffle the palette
#'   Default is `"normal"`.
#' @param palette Numeric or character specifying the palette. Options are
#'   1, 2, 3, 4, or palette name. Default is 1.
#' @param show_col Logical indicating whether to display the color palette.
#'   Default is `TRUE`.
#' @param seed Integer seed for random number generator when `cols = "random"`.
#'   Default is 123.
#'
#' @return Character vector of colors.
#'
#' @export
#'
#' @examples
#' # Get default palette
#' mycols <- get_cols()
#'
#' # Get random palette
#' mycols <- get_cols(cols = "random", seed = 456)
#'
#' # Use custom colors
#' mycols <- get_cols(cols = c("red", "blue", "green"))
get_cols <- function(cols = "normal",
                     palette = 1,
                     show_col = TRUE,
                     seed = 123) {
  # Map numeric palette to name
  palette_map <- c(
    "1" = "nrc",
    "2" = "jama",
    "3" = "aaas",
    "4" = "jco",
    "5" = "paired1",
    "6" = "paired2",
    "7" = "paired3",
    "8" = "paired4",
    "9" = "accent",
    "10" = "set2"
  )

  # Convert numeric palette to name
  if (is.numeric(palette) || grepl("^[0-9]+$", as.character(palette))) {
    palette_key <- as.character(as.numeric(palette))
    if (palette_key %in% names(palette_map)) {
      palette <- palette_map[palette_key]
    } else {
      palette <- "nrc" # default fallback
    }
  }

  # Handle custom color vector
  if (length(cols) > 1) {
    mycols <- cols
    if (show_col) {
      rlang::check_installed("scales")
      scales::show_col(mycols)
    }
    return(mycols)
  }

  # Handle single character specification
  cols <- as.character(cols)

  if (cols == "random") {
    # Get palette and shuffle
    mycols <- palettes(
      category = "random",
      palette = palette,
      show_col = FALSE,
      show_message = FALSE
    )

    cli::cli_alert_info("Using random seed: {seed}")
    set.seed(seed)
    mycols <- sample(mycols)

    if (show_col) {
      rlang::check_installed("scales")
      scales::show_col(mycols)
    }
  } else if (cols == "normal") {
    mycols <- palettes(
      category = "box",
      palette = palette,
      show_col = show_col,
      show_message = FALSE
    )
  } else {
    # Treat as palette name
    mycols <- palettes(
      category = "box",
      palette = cols,
      show_col = show_col,
      show_message = FALSE
    )
  }

  mycols
}
