#' Default Pattern List for Name Cleaning
#'
#' @description
#' A character vector of common substrings to remove from feature names.
#' Used in [remove_names()] and other helper functions.
#'
#' @format A character vector of length 12.
#'
#' @return Character vector of patterns to remove.
#'
#' @export
#'
#' @examples
#' # View default patterns
#' patterns_to_na
patterns_to_na <- c(
  "_cibersort", "xCell", "_EPIC", "_TIMER", "_quantiseq", "_MCP",
  "HALLMARK_", "_CIBERSORT", "xcell", "_timer", "_mcp", "_epic"
)


#' Remove Patterns from Column Names or Variables
#'
#' @description
#' Modifies column names or specified variables in a data frame by replacing
#' specified patterns with empty strings or spaces.
#'
#' @param input_df Data frame. Input data to modify.
#' @param variable Character. Column to modify: "colnames" for column names,
#'   or a specific column name. Default is "colnames".
#' @param patterns_to_na Character vector. Patterns to replace with empty
#'   string. Default uses [patterns_to_na].
#' @param patterns_space Character vector or `NULL`. Patterns to replace with
#'   spaces. Default is `NULL`.
#'
#' @return Modified data frame with patterns replaced.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' df <- data.frame(
#'   "CellA_cibersort" = 1:5,
#'   "CellB_xCell" = 6:10,
#'   "CellC_TIMER" = 11:15
#' )
#' result <- remove_names(df, variable = "colnames", patterns_to_na = patterns_to_na)
#' colnames(result)
remove_names <- function(input_df,
                         variable = "colnames",
                         patterns_to_na = patterns_to_na,
                         patterns_space = NULL) {
  if (!is.data.frame(input_df)) {
    cli::cli_abort("{.arg input_df} must be a data frame.")
  }

  if (variable != "colnames" && !variable %in% colnames(input_df)) {
    cli::cli_abort("Column {.val {variable}} not found in {.arg input_df}.")
  }

  if (variable == "colnames") {
    for (pattern in patterns_to_na) {
      colnames(input_df) <- gsub(colnames(input_df), pattern = pattern, replacement = "")
    }

    if (!is.null(patterns_space)) {
      for (pattern in patterns_space) {
        colnames(input_df) <- gsub(colnames(input_df), pattern = pattern, replacement = " ")
      }
    }
  } else {
    for (pattern in patterns_to_na) {
      input_df[[variable]] <- gsub(input_df[[variable]], pattern = pattern, replacement = "")
    }

    if (!is.null(patterns_space)) {
      for (pattern in patterns_space) {
        input_df[[variable]] <- gsub(input_df[[variable]], pattern = pattern, replacement = " ")
      }
    }
  }

  input_df
}
