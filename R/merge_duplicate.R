#' Merge Data Frames with Duplicated Column Names
#'
#' @description
#' Merges two data frames, resolving duplicated column names according to user
#' preference. Allows selection of which data frame's duplicated columns to
#' retain, ensuring data integrity during merging.
#'
#' @param x Data frame. First data frame to merge.
#' @param y Data frame. Second data frame to merge.
#' @param by.x Character. Column name(s) in `x` used for merging.
#' @param by.y Character. Column name(s) in `y` used for merging.
#' @param all.x Logical. Include all rows from `x` in output. Default is `FALSE`.
#' @param all.y Logical. Include all rows from `y` in output. Default is `FALSE`.
#' @param all Logical or `NULL`. If not `NULL`, include all rows from both `x`
#'   and `y`, overriding `all.x` and `all.y`.
#' @param choose Character. Which data frame's duplicated non-joining columns to
#'   retain: `"x"` or `"y"`. Default is `"x"`.
#'
#' @return Data frame resulting from merging `x` and `y` according to specified
#'   parameters.
#'
#' @export
#'
#' @examples
#' df1 <- data.frame(ID = 1:3, Name = c("A", "B", "C"), Value = 1:3)
#' df2 <- data.frame(ID = 1:3, Name = c("X", "Y", "Z"), Score = 4:6)
#'
#' # Merge keeping duplicated columns from x
#' merged_df <- merge_duplicate(df1, df2,
#'   by.x = "ID", by.y = "ID",
#'   all.x = TRUE, choose = "x"
#' )
#' print(merged_df)
#'
#' # Merge keeping duplicated columns from y
#' merged_df2 <- merge_duplicate(df1, df2,
#'   by.x = "ID", by.y = "ID",
#'   all = TRUE, choose = "y"
#' )
merge_duplicate <- function(x, y, by.x, by.y,
                            all.x = FALSE, all.y = FALSE,
                            all = NULL,
                            choose = c("x", "y")) {
  # Input validation
  if (!is.data.frame(x)) {
    cli::cli_abort("{.arg x} must be a data frame.")
  }
  if (!is.data.frame(y)) {
    cli::cli_abort("{.arg y} must be a data frame.")
  }

  choose <- rlang::arg_match(choose)

  if (!by.x %in% colnames(x)) {
    cli::cli_abort("Column {.val {by.x}} not found in {.arg x}.")
  }
  if (!by.y %in% colnames(y)) {
    cli::cli_abort("Column {.val {by.y}} not found in {.arg y}.")
  }

  # Find duplicate column names (excluding merge keys)
  duplicate_names <- intersect(colnames(x), colnames(y))

  if (by.x == by.y) {
    duplicate_names <- setdiff(duplicate_names, by.x)
  } else {
    duplicate_names <- setdiff(duplicate_names, c(by.x, by.y))
  }

  # Remove duplicates from the non-chosen data frame
  if (length(duplicate_names) > 0) {
    cli::cli_alert_info(
      "Removing {length(duplicate_names)} duplicate column{?s} from {.arg {ifelse(choose == 'x', 'y', 'x')}}"
    )

    if (choose == "x") {
      y <- y[, !colnames(y) %in% duplicate_names, drop = FALSE]
    } else {
      x <- x[, !colnames(x) %in% duplicate_names, drop = FALSE]
    }
  }

  # Perform merge
  if (!is.null(all)) {
    res <- merge(x = x, y = y, by.x = by.x, by.y = by.y, all = all)
  } else {
    res <- merge(x = x, y = y, by.x = by.x, by.y = by.y, all.x = all.x, all.y = all.y)
  }

  res
}
