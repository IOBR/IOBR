#' Create Nested Output Folders
#'
#' Creates one to three nested folders (if not existing) under the current working directory and returns their names and absolute path.
#'
#' @param f1 Character. First-level folder name.
#' @param f2 Character or NULL. Second-level folder name.
#' @param f3 Character or NULL. Third-level folder name.
#' @param return Deprecated (not used). Kept for backward compatibility.
#'
#' @return List with elements: `folder_name` (relative) and `abspath` (absolute path ending with '/').
#' @export
#' @examples
#' creat_folder("1-result")
creat_folder <- function(f1, f2 = NULL, f3 = NULL, return = NULL) {
  if (!is.null(f3)) {
    path <- file.path(getwd(), f1, f2, f3)
    if (!dir.exists(path)) dir.create(file.path(getwd(), f1, f2, f3), recursive = TRUE)
  } else if (!is.null(f2)) {
    path <- file.path(getwd(), f1, f2)
    if (!dir.exists(path)) dir.create(file.path(getwd(), f1, f2), recursive = TRUE)
  } else {
    path <- file.path(getwd(), f1)
    if (!dir.exists(path)) dir.create(file.path(getwd(), f1), recursive = TRUE)
  }


  if (is.null(return)) {
    if (!is.null(f3)) {
      res <- list(
        "folder_name" = paste0(f1, "/", f2, "/", f3),
        "abspath" = paste0(file.path(getwd(), f1, f2, f3), "/")
      )
    } else if (!is.null(f2)) {
      res <- list(
        "folder_name" = paste0(f1, "/", f2),
        "abspath" = paste0(file.path(getwd(), f1, f2), "/")
      )
    } else {
      res <- list(
        "folder_name" = f1,
        "abspath" = paste0(file.path(getwd(), f1), "/")
      )
    }
  } else {
    if (return == 1) {
      res <- list(
        "folder_name" = f1,
        "abspath" = paste0(file.path(getwd(), f1), "/")
      )
    } else if (return == 2) {
      res <- list(
        "folder_name" = paste0(f1, "/", f2),
        "abspath" = paste0(file.path(getwd(), f1, f2), "/")
      )
    } else if (return == 3) {
      res <- list(
        "folder_name" = paste0(f1, "/", f2, "/", f3),
        "abspath" = paste0(file.path(getwd(), f1, f2, f3), "/")
      )
    }
  }

  return(res)
}
