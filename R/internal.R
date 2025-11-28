# Internal helper to load package data safely
.load_data <- function(name) {
  utils::data(list = name, package = "IOBR", envir = environment())
  get(name, envir = environment(), inherits = FALSE)
}
