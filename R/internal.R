#' GSVA API Version Detection
#'
#' @description
#' Detects whether the installed GSVA package supports the new
#' parameter-based API (gsvaParam/ssgseaParam) or the old direct argument API.
#' This function is used internally to ensure compatibility across different
#' GSVA package versions (1.50.0+ vs older versions).
#'
#' @return A list with two elements:
#'   - `use_new_api`: Logical indicating whether to use the new API (`TRUE`)
#'     or old API (`FALSE`)
#'   - `gsva_version`: Character string of the installed GSVA version,
#'     or `"not installed"` if not available
#'
#' @keywords internal
gsva_use_new_api <- function() {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    return(list(use_new_api = FALSE, gsva_version = "not installed"))
  }

  gsva_version <- as.character(utils::packageVersion("GSVA"))

  # Check for new API functions (introduced in GSVA 1.50.0)
  use_new_api <- exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE) &&
    exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)

  list(use_new_api = use_new_api, gsva_version = gsva_version)
}

#' Load IOBR Datasets
#'
#' @description
#' Loads internal datasets from the IOBR package. Supports both sysdata (internal)
#' and exported data files included in the package.
#'
#' @param name Character string. Name of the dataset to load. Must be a single value.
#'   Available datasets include:
#'   - Expression data: `"eset_stad"`, `"imvigor210_eset"`, `"melanoma_data"`
#'   - Signatures: `"signature_tme"`, `"signature_metabolism"`, `"signature_collection"`
#'   - Gene sets: `"hallmark"`, `"kegg"`, `"go_bp"`, `"go_cc"`, `"go_mf"`
#'   - Cell markers: `"cellmarkers"`, `"mcp_genes"`
#'   - Phenotype data: `"pdata_stad"`, `"pdata_sig_tme"`, `"pdata_acrg"`
#'   - Reference data: `"xCell.data"`, `"quantiseq_data"`, `"TRef"`, `"BRef"`
#'   - Color palettes: `"palette1"`, `"palette2"`, `"palette3"`, `"palette4"`
#'
#' @returns Dataset object, typically a `list`, `data.frame`, or `matrix`.
#'   The exact type depends on the requested dataset.
#'
#' @export
#'
#' @examples
#' # Load signature collection (stored in sysdata, no download)
#' sig_tme <- load_data("signature_tme")
#'
#' # Load color palette (stored in sysdata, no download)
#' colors <- load_data("palette1")
#'
#' # Error handling with suggestions for similar names
#' try(load_data("sign_tme")) # Will suggest "signature_tme"
#'
#' if (interactive()) {
#'   # Load expression data (triggers download from GitHub)
#'   eset <- load_data("eset_stad")
#' }
load_data <- function(name) {
  # Input validation
  if (!is.character(name) || length(name) != 1 || is.na(name) || nchar(name) == 0) {
    stop("'name' must be a single non-empty character string")
  }

  # Define sysdata names (internal package data - minimal essential data)
  # Large datasets moved to GitHub to meet CRAN size requirements
  sysdata_names <- c(
    # Core signatures (small, frequently used) - stored in sysdata.rda
    "signature_tme",

    # MCPCounter (required for deconvolution) - stored in sysdata.rda
    "mcp_genes", "mcp_probesets",

    # Color palettes and configs (small) - stored in sysdata.rda
    "palette1", "palette2", "palette3", "palette4",
    "panel_for_gene", "panel_for_signature",
    "patterns_to_na", "sig_excel"
  )

  # Get available data names from package
  pkg_data <- utils::data(package = "IOBR")$results
  data_names <- if (nrow(pkg_data) > 0) pkg_data[, "Item"] else character(0)

  # Load data based on type
  if (name %in% sysdata_names) {
    # Internal sysdata - use eval for lazy loading
    result <- eval(parse(text = name))
    return(result)
  }

  if (name %in% data_names) {
    # External data file
    data(list = name, package = "IOBR", envir = environment())
    return(get(name, envir = environment()))
  }

  # Check if it's a GitHub-hosted dataset (large datasets moved from sysdata and data/)
  github_datasets <- list_github_datasets()

  if (name %in% github_datasets) {
    # Try to download from GitHub
    return(download_iobr_data(name))
  }

  # Dataset not found - provide helpful error message with suggestions
  available <- sort(c(sysdata_names, data_names, github_datasets))

  # Find similar names for suggestions
  distances <- utils::adist(name, available, ignore.case = TRUE)
  suggestions <- available[which(distances <= 3)]

  error_msg <- c(
    "Dataset {.val {name}} not found in IOBR package.",
    "i" = paste(
      "Available datasets:",
      paste(utils::head(available, 100), collapse = ", "),
      ifelse(length(available) > 100, "...", "")
    )
  )

  if (length(suggestions) > 0) {
    error_msg <- c(
      error_msg,
      "!" = "Did you mean: {.val {suggestions}}?"
    )
  }

  cli::cli_abort(error_msg)
}
