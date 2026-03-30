#' MCP-counter Cell Population Abundance Estimation
#'
#' @description
#' Estimates the abundance of different immune and stromal cell populations
#' using the MCP-counter method. Works with various gene identifiers including
#' Affymetrix probesets, HUGO gene symbols, Entrez IDs, and Ensembl IDs.
#'
#' @param expression Matrix or data.frame with features in rows and samples
#'   in columns.
#' @param featuresType Type of identifiers for expression features.
#'   Options: "affy133P2_probesets", "HUGO_symbols", "ENTREZ_ID", "ENSEMBL_ID".
#'   Default is "affy133P2_probesets".
#' @param probesets Probesets data table. Default loads from GitHub.
#' @param genes Genes data table. Default loads from GitHub.
#'
#' @return Matrix with cell populations in rows and samples in columns.
#'
#' @author Etienne Becht
#' @export
#'
#' @examples
#' \donttest{
#' # Create example expression data
#' expr <- matrix(runif(1000), nrow = 100, ncol = 10)
#' rownames(expr) <- paste0("Gene", 1:100)
#' # Estimate with HUGO symbols
#' estimates <- MCPcounter.estimate(expr, featuresType = "HUGO_symbols")
#' }
MCPcounter.estimate <- function(
    expression,
    featuresType = c("affy133P2_probesets", "HUGO_symbols", "ENTREZ_ID", "ENSEMBL_ID"),
    probesets = read.table(
      url(paste0(
        "https://raw.githubusercontent.com/ebecht/",
        "MCPcounter/master/Signatures/probesets.txt"
      )),
      sep = "\t", stringsAsFactors = FALSE, colClasses = "character"
    ),
    genes = read.table(
      url(paste0(
        "https://raw.githubusercontent.com/ebecht/",
        "MCPcounter/master/Signatures/genes.txt"
      )),
      sep = "\t", stringsAsFactors = FALSE, header = TRUE,
      colClasses = "character", check.names = FALSE
    )
) {
  featuresType <- rlang::arg_match(featuresType)

  # Get features based on type
  result <- switch(featuresType,
    affy133P2_probesets = .get_probeset_features(expression, probesets),
    HUGO_symbols = .get_hugo_features(expression, genes),
    ENTREZ_ID = .get_entrez_features(expression, genes),
    ENSEMBL_ID = .get_ensembl_features(expression, genes)
  )

  features <- result$features
  missing.populations <- result$missing

  if (length(missing.populations) > 0) {
    cli::cli_warn("No markers for population(s): {.val {missing.populations}}")
  }

  t(.appendSignatures(expression, features))
}

#' @keywords internal
.get_probeset_features <- function(expression, probesets) {
  markers.names <- unique(probesets[, 2])
  features <- split(probesets[, 1], probesets[, 2])
  features <- lapply(features, intersect, x = rownames(expression))
  features <- features[vapply(features, length, integer(1)) > 0]
  missing.populations <- setdiff(markers.names, names(features))
  features <- features[intersect(markers.names, names(features))]

  list(features = features, missing = missing.populations)
}

#' @keywords internal
.get_hugo_features <- function(expression, genes) {
  markersG <- genes
  features <- subset(markersG, get("HUGO symbols") %in% rownames(expression))
  markers.names <- unique(features[, "Cell population"])
  features <- split(features[, "HUGO symbols"], features[, "Cell population"])
  missing.populations <- setdiff(markers.names, names(features))
  features <- features[intersect(markers.names, names(features))]

  list(features = features, missing = missing.populations)
}

#' @keywords internal
.get_entrez_features <- function(expression, genes) {
  markersG <- genes
  features <- subset(markersG, ENTREZID %in% rownames(expression))
  markers.names <- unique(features[, "Cell population"])
  features <- split(features[, "ENTREZID"], features[, "Cell population"])
  missing.populations <- setdiff(markers.names, names(features))
  features <- features[intersect(markers.names, names(features))]

  list(features = features, missing = missing.populations)
}

#' @keywords internal
.get_ensembl_features <- function(expression, genes) {
  markersG <- genes
  features <- subset(markersG, get("ENSEMBL ID") %in% rownames(expression))
  markers.names <- unique(features[, "Cell population"])
  features <- split(features[, "ENSEMBL ID"], features[, "Cell population"])
  missing.populations <- setdiff(markers.names, names(features))
  features <- features[intersect(markers.names, names(features))]

  list(features = features, missing = missing.populations)
}

#' Append Signatures to Expression Matrix
#'
#' @description
#' Calculates mean expression for each marker feature set.
#'
#' @param xp Expression matrix with features in rows and samples in columns.
#' @param markers List of marker gene vectors for each cell population.
#'
#' @return Matrix with summarized expression values.
#'
#' @keywords internal
.appendSignatures <- function(xp, markers) {
  res <- as.data.frame(do.call(
    cbind,
    lapply(markers, function(x) {
      common_genes <- intersect(row.names(xp), x)
      if (length(common_genes) == 0) {
        return(rep(NA, ncol(xp)))
      }
      apply(xp[common_genes, , drop = FALSE], 2, mean, na.rm = TRUE)
    })
  ))
  res
}

#' Test for Cell Population Infiltration
#'
#' @description
#' Returns p-values for the null hypothesis that samples are not infiltrated
#' by the corresponding cell population.
#'
#' @param MCPcounterMatrix Matrix, usually output from MCPcounter.estimate.
#' @param platform Expression platform: "133P2", "133A", or "HG1".
#'   Default is "133P2".
#'
#' @return Matrix with samples in rows and cell populations in columns.
#'   Elements are p-values.
#'
#' @author Etienne Becht
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data
#' scores <- matrix(runif(30), nrow = 3, ncol = 10)
#' rownames(scores) <- c("T cells", "B cells", "NK cells")
#' pvals <- test_for_infiltration(scores, platform = "133P2")
#' }
test_for_infiltration <- function(MCPcounterMatrix,
                                  platform = c("133P2", "133A", "HG1")) {
  platform <- rlang::arg_match(platform)

  MCPcounterMatrix <- t(MCPcounterMatrix)
  params <- null_models[grep(platform, colnames(null_models))]
  rownames(params) <- null_models[, "Cell.population"]
  colnames(params) <- sub(platform, "", colnames(params), fixed = TRUE)

  res <- vapply(colnames(MCPcounterMatrix), function(x) {
    stats::pnorm(
      MCPcounterMatrix[, x],
      mean = params[x, "mu."],
      sd = params[x, "sigma."],
      lower.tail = FALSE
    )
  }, numeric(nrow(MCPcounterMatrix)))

  rownames(res) <- rownames(MCPcounterMatrix)
  res
}
