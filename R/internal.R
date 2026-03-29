.iobr_cache <- new.env(parent = emptyenv())

.load_data <- function(name) {
  cache <- .iobr_cache

  # 1. cache
  if (exists(name, envir = cache, inherits = FALSE)) {
    return(get(name, envir = cache, inherits = FALSE))
  }

  # 2. namespace (for sysdata.rda or internal objects)
  ns <- asNamespace("IOBR")
  if (exists(name, envir = ns, inherits = FALSE)) {
    obj <- get(name, envir = ns, inherits = FALSE)
    assign(name, obj, envir = cache)
    return(obj)
  }

  # 3. data/
  tmp <- new.env(parent = emptyenv())
  tryCatch(
    {
      utils::data(list = name, package = "IOBR", envir = tmp)

      if (!exists(name, envir = tmp, inherits = FALSE)) {
        stop("not found")
      }

      obj <- get(name, envir = tmp, inherits = FALSE)
      assign(name, obj, envir = cache)
      obj
    },
    error = function(e) {
      available_cache <- ls(envir = cache, all.names = TRUE)
      available_cache <- available_cache[!startsWith(available_cache, ".")]

      available_data <- tryCatch({
        info <- utils::data(package = "IOBR")
        if (!is.null(info$results) && "Item" %in% colnames(info$results)) {
          unique(info$results[, "Item"])
        } else {
          character(0)
        }
      }, error = function(e) character(0))

      available_ns <- ls(envir = ns, all.names = TRUE)
      available_ns <- available_ns[!startsWith(available_ns, ".")]

      available <- sort(unique(c(available_cache, available_ns, available_data)))

      stop(
        "Data '", name, "' not found in IOBR. Available objects/datasets: ",
        paste(available, collapse = ", "),
        call. = FALSE
      )
    }
  )
}