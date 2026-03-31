#' Create Nested Output Folders
#'
#' @description
#' Creates one to three nested folders (if not existing) under the current
#' working directory and returns their names and absolute paths.
#'
#' @param f1 Character. First-level folder name.
#' @param f2 Character or `NULL`. Second-level folder name. Default is `NULL`.
#' @param f3 Character or `NULL`. Third-level folder name. Default is `NULL`.
#' @param return Deprecated (not used). Kept for backward compatibility.
#'
#' @return List with elements:
#' \describe{
#'   \item{folder_name}{Relative path to the created folder}
#'   \item{abspath}{Absolute path ending with '/'}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Create single folder (in temp directory for examples)
#' oldwd <- getwd()
#' on.exit(setwd(oldwd))
#' setwd(tempdir())
#' creat_folder("1-result")
#'
#' # Create nested folders
#' creat_folder("1-result", "figures", "correlation")
#' }
creat_folder <- function(f1, f2 = NULL, f3 = NULL, return = NULL) {
  # Validate inputs
  if (!is.character(f1) || length(f1) != 1 || nchar(f1) == 0) {
    cli::cli_abort("{.arg f1} must be a non-empty character string")
  }

  # Build path components
  components <- c(f1, f2, f3)
  components <- components[!vapply(components, is.null, logical(1))]

  # Create full path
  path <- do.call(file.path, c(list(getwd()), components))

  # Create directory if it doesn't exist
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  # Build return values
  folder_name <- paste(components, collapse = "/")
  abspath <- paste0(path, "/")

  # Handle deprecated 'return' parameter
  if (!is.null(return)) {
    cli::cli_warn("{.arg return} parameter is deprecated and will be ignored")
  }

  list(
    folder_name = folder_name,
    abspath = abspath
  )
}
