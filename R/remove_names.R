#' @export
patterns_to_na <- c("_cibersort", "xCell", "_EPIC", "_TIMER", "_quantiseq", "_MCP", "HALLMARK_", "_CIBERSORT", "xcell", "_timer", "_mcp", "_epic")



#' Remove Patterns from Column Names or Variables
#'
#' This function modifies column names or specified variables in a data frame by replacing
#' specified patterns with NA or spaces.
#'
#' @param input_df Input data frame.
#' @param variable Column to modify: "colnames" for column names, or a specific column name. Default is "colnames".
#' @param patterns_to_na Vector of patterns to replace with empty string. Default uses predefined patterns.
#' @param patterns_space Vector of patterns to replace with spaces. Default is NULL.
#'
#' @return Modified data frame with patterns replaced.
#' @export
#' @examples
#' data("imvigor210_sig", package = "IOBR")
#' input <- remove_names(imvigor210_sig, variable = "colnames", patterns_to_na = patterns_to_na, patterns_space = NULL)
remove_names <- function(input_df, variable = "colnames", patterns_to_na = patterns_to_na, patterns_space = NULL) {
  if (variable == "colnames") {
    for (i in 1:length(patterns_to_na)) {
      colnames(input_df) <- gsub(colnames(input_df), pattern = patterns_to_na[i], replacement = "")
    }

    if (!is.null(patterns_space)) {
      for (j in 1:length(patterns_space)) {
        colnames(input_df) <- gsub(colnames(input_df), pattern = patterns_space[j], replacement = " ")
      }
    }
  } else {
    input_df <- as.data.frame(input_df)
    for (i in 1:length(patterns_to_na)) {
      input_df[, variable] <- gsub(input_df[, variable], pattern = patterns_to_na[i], replacement = "")
    }
    if (!is.null(patterns_space)) {
      for (j in 1:length(patterns_space)) {
        input_df[, variable] <- gsub(input_df[, variable], pattern = patterns_space[j], replacement = " ")
      }
    }
  }
  return(input_df)
}
