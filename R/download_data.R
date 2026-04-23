#' Download IOBR Data from GitHub with Mirror Support
#'
#' @description
#' Downloads large datasets from GitHub releases to avoid CRAN size limits.
#' Supports multiple download mirrors for users in different regions.
#' Data is cached locally after first download. Cache directory can be
#' customized via `options(IOBR.cache_dir = "your/path")`.
#'
#' @param name Character string. Name of the dataset to download.
#' @param force Logical. Whether to force re-download even if cached. Default: FALSE.
#' @param verbose Logical. Whether to print progress messages. Default: TRUE.
#' @param mirrors Character vector. URLs of mirrors to try. Default uses
#'   get_default_mirrors().
#' @param cache_dir Character string. Custom cache directory. If NULL, uses
#'   the option `IOBR.cache_dir` or the default system cache location.
#'
#' @return The requested dataset.
#' @export
#'
#' @examples
#' \donttest{
#' # Download TCGA STAD signature data
#' tcga_sig <- download_iobr_data("tcga_stad_sig")
#'
#' # Download with custom mirrors
#' eset <- download_iobr_data("eset_stad",
#'   mirrors = c(
#'     "https://ghproxy.vip/https://github.com",
#'     "https://gh-proxy.org/https://github.com"
#'   )
#' )
#'
#' # Use custom cache directory (use tempdir() for examples)
#' options(IOBR.cache_dir = tempdir())
#' data <- download_iobr_data("lm22")
#' }
download_iobr_data <- function(name, force = FALSE, verbose = TRUE,
                               mirrors = get_default_mirrors(),
                               cache_dir = NULL) {
  # Get all available GitHub datasets
  github_data <- list_github_datasets()

  if (!name %in% github_data) {
    stop(sprintf(
      "Dataset '%s' not available for download. Available: %s",
      name, paste(github_data, collapse = ", ")
    ))
  }

  # Set up cache directory (priority: argument > option > default)
  cache_dir <- get_iobr_cache_dir(cache_dir)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cache_file <- file.path(cache_dir, sprintf("%s.rda", name))

  # Check if already cached
  if (file.exists(cache_file) && !force) {
    if (verbose) cli::cli_alert_info("Loading cached data: {.val {name}}")
    env <- new.env()
    load(cache_file, envir = env)
    # Return the object (might be named 'data' or the actual name)
    obj_names <- ls(env)
    if (length(obj_names) == 1) {
      return(env[[obj_names[1]]])
    } else {
      return(env[[name]])
    }
  }

  # Try each mirror until success
  last_error <- NULL
  for (i in seq_along(mirrors)) {
    mirror <- mirrors[i]
    url <- sprintf(
      "%s/IOBR/IOBR/releases/download/data-v1.0/%s.rda",
      mirror, name
    )

    if (verbose) {
      cli::cli_alert_info("Trying mirror {i}/{length(mirrors)}: {.url {mirror}}")
    }

    tryCatch(
      {
        utils::download.file(url, cache_file, mode = "wb", quiet = !verbose)

        if (file.exists(cache_file) && file.size(cache_file) > 0) {
          if (verbose) cli::cli_alert_success("Download complete: {.val {name}}")

          env <- new.env()
          load(cache_file, envir = env)
          obj_names <- ls(env)
          if (length(obj_names) == 1) {
            return(env[[obj_names[1]]])
          } else {
            return(env[[name]])
          }
        }
      },
      error = function(e) {
        last_error <<- e
        if (verbose) cli::cli_alert_warning("Mirror {i} failed: {e$message}")
        # Clean up partial download
        if (file.exists(cache_file)) file.remove(cache_file)
      }
    )
  }

  # All mirrors failed
  cache_dir <- get_iobr_cache_dir()
  manual_url <- sprintf(
    "https://github.com/IOBR/IOBR/releases/download/data-v1.0/%s.rda", name
  )
  cli::cli_alert_danger("All download mirrors failed for dataset: {.val {name}}")
  cli::cli_alert_info("Please try the following manual download steps:")
  cli::cli_ul(c(
    "1. Download the file from: {.url {manual_url}}",
    sprintf("2. Save it to: {.path %s}", file.path(cache_dir, paste0(name, ".rda"))),
    "3. Run your code again - the data will be loaded from cache"
  ))
  stop(sprintf(
    "Failed to download '%s' from all %d mirrors. Please download manually.",
    name, length(mirrors)
  ))
}

#' Get Default Download Mirrors
#'
#' @description Returns the default list of download mirrors.
#' @return Character vector of mirror URLs.
#' @keywords internal
get_default_mirrors <- function() {
  c(
    # Original GitHub (default)
    "https://github.com",
    # Chinese mirrors
    "https://ghproxy.vip/https://github.com",
    "https://gh-proxy.org/https://github.com",
    "https://ghfast.top/https://github.com"
  )
}

#' Add Custom Download Mirror
#'
#' @description
#' Adds a custom mirror URL to the default mirrors for the current session.
#' The mirror URL should be a base URL that will be prepended to GitHub paths.
#'
#' @param url Character string. The mirror URL to add.
#' @param position Character. Where to add the mirror: "first", "last", or
#'   "before_github". Default: "first".
#'
#' @return Invisibly returns the updated mirror list.
#' @export
#'
#' @examples
#' \donttest{
#' # Add a custom mirror to try first
#' add_iobr_mirror("https://my-mirror.com/https://github.com")
#'
#' # Add mirror to try before default GitHub
#' add_iobr_mirror("https://fast-mirror.org", position = "before_github")
#'
#' # Download with the new mirror
#' data <- download_iobr_data("BRef")
#' }
add_iobr_mirror <- function(url, position = c("first", "last", "before_github")) {
  position <- match.arg(position)

  # Validate URL
  if (!grepl("^https?://", url)) {
    stop("Invalid URL. Must start with http:// or https://")
  }

  # Get current mirrors from option or default
  current_mirrors <- getOption("IOBR.download_mirrors", get_default_mirrors())

  # Remove trailing slash if present
  url <- sub("/$", "", url)

  # Add URL if not already present
  if (url %in% current_mirrors) {
    cli::cli_alert_info("Mirror {.url {url}} already exists in the list")
    return(invisible(current_mirrors))
  }

  # Add at specified position
  new_mirrors <- switch(position,
    "first" = c(url, current_mirrors),
    "last" = c(current_mirrors, url),
    "before_github" = {
      # Insert before the default GitHub URL
      github_idx <- which(current_mirrors == "https://github.com")
      if (length(github_idx) == 0) {
        c(current_mirrors, url)
      } else {
        c(current_mirrors[1:(github_idx - 1)], url, current_mirrors[github_idx:length(current_mirrors)])
      }
    }
  )

  # Store in options
  options(IOBR.download_mirrors = new_mirrors)

  cli::cli_alert_success("Added mirror {.url {url}} to position: {.val {position}}")
  cli::cli_alert_info("Current mirrors: {.val {length(new_mirrors)}} total")

  invisible(new_mirrors)
}

#' List Current Download Mirrors
#'
#' @description Returns the current list of download mirrors.
#' @return Character vector of mirror URLs.
#' @export
#'
#' @examples
#' list_iobr_mirrors()
list_iobr_mirrors <- function() {
  mirrors <- getOption("IOBR.download_mirrors", get_default_mirrors())

  cat("Current IOBR download mirrors:\n")
  cat("==============================\n")
  for (i in seq_along(mirrors)) {
    cat(sprintf("%d. %s\n", i, mirrors[i]))
  }
  cat("==============================\n")
  cat(sprintf("Total: %d mirrors\n", length(mirrors)))

  invisible(mirrors)
}

#' Reset Download Mirrors to Default
#'
#' @description Resets the download mirror list to the default values.
#' @return Invisibly returns the default mirror list.
#' @export
#'
#' @examples
#' reset_iobr_mirrors()
reset_iobr_mirrors <- function() {
  options(IOBR.download_mirrors = NULL)
  cli::cli_alert_success("Download mirrors reset to default")
  invisible(get_default_mirrors())
}

#' Get IOBR Cache Directory
#'
#' @description Returns the current cache directory for IOBR downloaded data.
#' The cache directory is determined in the following priority order:
#' 1. Function argument `cache_dir` (if provided)
#' 2. Option `IOBR.cache_dir` (if set via `options()`)
#' 3. Default system cache location via `tools::R_user_dir()`
#'
#' @param cache_dir Optional character string to override the current setting.
#' @return Character string with the cache directory path.
#' @export
#'
#' @examples
#' # Get current cache directory
#' get_iobr_cache_dir()
#'
#' # Set custom cache directory via options (use tempdir() for examples)
#' options(IOBR.cache_dir = tempdir())
#' get_iobr_cache_dir()
get_iobr_cache_dir <- function(cache_dir = NULL) {
  # Priority: argument > option > default
  if (!is.null(cache_dir)) {
    if (!is.character(cache_dir) || length(cache_dir) != 1) {
      stop("'cache_dir' must be a single character string")
    }
    return(path.expand(cache_dir))
  }

  option_cache <- getOption("IOBR.cache_dir")
  if (!is.null(option_cache)) {
    return(path.expand(option_cache))
  }

  # Default to R user cache directory
  tools::R_user_dir("IOBR", which = "cache")
}

#' Set IOBR Cache Directory
#'
#' @description Sets a custom cache directory for IOBR downloaded data.
#' This is useful when you want to store cached data in a specific location,
#' such as a shared network drive or a custom directory.
#'
#' @param path Character string. The path to the cache directory.
#' @param create Logical. Whether to create the directory if it doesn't exist.
#'   Default: TRUE.
#'
#' @return Invisibly returns the cache directory path.
#' @export
#'
#' @examples
#' \donttest{
#' # Set a custom cache directory (use tempdir() for examples)
#' set_iobr_cache_dir(tempdir())
#'
#' # Check the current cache directory
#' get_iobr_cache_dir()
#'
#' # Download data will now use the custom cache
#' data <- download_iobr_data("lm22")
#' }
set_iobr_cache_dir <- function(path, create = TRUE) {
  if (!is.character(path) || length(path) != 1) {
    stop("'path' must be a single character string")
  }

  path <- path.expand(path)

  if (create && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) {
      stop(sprintf("Failed to create cache directory: %s", path))
    }
    cli::cli_alert_success("Created cache directory: {.path {path}}")
  }

  options(IOBR.cache_dir = path)
  cli::cli_alert_success("IOBR cache directory set to: {.path {path}}")

  invisible(path)
}

#' Reset IOBR Cache Directory to Default
#'
#' @description Resets the cache directory to the default system location.
#' @return Invisibly returns the default cache directory path.
#' @export
#'
#' @examples
#' reset_iobr_cache_dir()
reset_iobr_cache_dir <- function() {
  options(IOBR.cache_dir = NULL)
  default_cache <- tools::R_user_dir("IOBR", which = "cache")
  cli::cli_alert_success("Cache directory reset to default: {.path {default_cache}}")
  invisible(default_cache)
}

#' List Available GitHub Datasets
#'
#' @return Character vector of available dataset names.
#' @export
#'
#' @examples
#' list_github_datasets()
list_github_datasets <- function() {
  c(
    # From data/ directory - Reference matrices
    "BRef", "TRef", "lm22",

    # From data/ directory - Annotation files
    "anno_gc_vm32", "anno_grch38", "anno_hug133plus2", "anno_illumina", "anno_rnaseq",

    # From data/ directory - Example datasets
    "tcga_stad_sig", "imvigor210_sig", "eset_stad", "sig_stad",
    "eset_gse62254", "eset_tme_stad", "eset_blca", "deg",

    # From sysdata - Gene sets
    "hallmark", "kegg", "go_bp", "go_cc", "go_mf", "reactome",
    "msig_immune", "msig_sc",

    # From sysdata - Reference data
    "xCell.data", "quantiseq_data",

    # From sysdata - Cell markers and gene lists
    "cancer_type_genes", "cellmarkers", "common_genes", "immuneCuratedData",
    "ips_gene_set", "length_ensembl", "mRNA_cell_default",
    "mus_human_gene_symbol", "onco_sig", "PurityDataAffy", "SI_geneset",

    # From sysdata - Signatures
    "signature_collection_citation", "signature_metabolism",
    "signature_sc", "signature_tumor",

    # From sysdata - Example datasets
    "imvigor210_eset", "melanoma_data",
    "pdata_acrg", "pdata_GSE63557", "pdata_sig_tme",
    "pdata_sig_tme_binary", "pdata_tme_binary", "tcga_stad_var"
  )
}

#' Clear IOBR Data Cache
#'
#' @description Removes all cached data files downloaded from GitHub.
#' @param cache_dir Character string. Custom cache directory. If NULL, uses
#'   the option `IOBR.cache_dir` or the default system cache location.
#' @return Invisible NULL. Called for side effects of clearing the cache.
#' @export
#'
#' @examples
#' clear_iobr_cache()
clear_iobr_cache <- function(cache_dir = NULL) {
  cache_dir <- get_iobr_cache_dir(cache_dir)
  if (dir.exists(cache_dir)) {
    files <- list.files(cache_dir, full.names = TRUE)
    if (length(files) > 0) {
      file.remove(files)
      cli::cli_alert_success("Cache cleared: {.val {length(files)}} file(s) removed from {.path {cache_dir}}")
    } else {
      cli::cli_alert_info("Cache is already empty: {.path {cache_dir}}")
    }
  } else {
    cli::cli_alert_info("Cache directory does not exist: {.path {cache_dir}}")
  }
  invisible(NULL)
}
