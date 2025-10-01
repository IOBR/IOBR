#' Append Signatures to Expression Matrix
#'
#' Takes as input an expression matrix and a list of marker features and returns summarized expression values.
#'
#' @param xp An expression matrix with features in rows and samples in columns.
#' @param markers A list whose names are cellular populations' names and elements are character vectors of features.
#'
#' @return Matrix with the summarized expression of each marker feature set in rows.
#' @author Etienne Becht
appendSignatures <- function(xp, markers) {
  res <- as.data.frame(do.call(
    cbind,
    lapply(markers, function(x) {
      apply(xp[intersect(row.names(xp), x), , drop = F], 2, mean, na.rm = T)
    })
  ))
  res
}

#' MCP-counter Cell Population Abundance Estimation
#'
#' Produces a matrix with abundance estimates from an expression matrix using MCP-counter method.
#'
#' @param expression Matrix or data.frame with features in rows and samples in columns.
#' @param featuresType Type of identifiers for expression features. Defaults to "affy133P2_probesets" for Affymetrix Human Genome 133 Plus 2.0 probesets. Other options are "HUGO_symbols" (Official gene symbols), "ENTREZ_ID" (Entrez Gene ID) or "ENSEMBL_ID" (ENSEMBL Gene ID).
#' @param probesets Probesets data table (default loads from GitHub).
#' @param genes Genes data table (default loads from GitHub).
#'
#' @return Matrix with cell populations in rows and samples in columns.
#' @author Etienne Becht
#' @examples
#' # Example usage (requires appropriate expression data)
#' # estimates <- MCPcounter.estimate(expression_matrix, featuresType = "HUGO_symbols")
MCPcounter.estimate <- function(
    expression,
    featuresType = c("affy133P2_probesets", "HUGO_symbols", "ENTREZ_ID", "ENSEMBL_ID")[1],
    probesets = read.table(curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"), sep = "\t", stringsAsFactors = FALSE, colClasses = "character"),
    genes = read.table(curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"), sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)) {
  ## marker.names=c("T cells","CD8 T cells","Cytotoxic lymphocytes","NK cells","B lineage","Monocytic lineage","Myeloid dendritic cells","Neutrophils","Endothelial cells","Fibroblasts")


  if (featuresType == "affy133P2_probesets") {
    features <- probesets
    markers.names <- unique(features[, 2])
    features <- split(features[, 1], features[, 2])
    features <- lapply(features, intersect, x = rownames(expression))
    features <- features[sapply(features, function(x) length(x) > 0)]
    missing.populations <- setdiff(markers.names, names(features))
    features <- features[intersect(markers.names, names(features))]
  } else {
    markersG <- genes
  }

  if (featuresType == "HUGO_symbols") {
    features <- subset(markersG, get("HUGO symbols") %in% rownames(expression))
    markers.names <- unique(features[, "Cell population"])
    features <- split(features[, "HUGO symbols"], features[, "Cell population"])
    missing.populations <- setdiff(markers.names, names(features))
    features <- features[intersect(markers.names, names(features))]
  }

  if (featuresType == "ENTREZ_ID") {
    features <- subset(markersG, ENTREZID %in% rownames(expression))
    markers.names <- unique(features[, "Cell population"])
    features <- split(features[, "ENTREZID"], features[, "Cell population"])
    missing.populations <- setdiff(markers.names, names(features))
    features <- features[intersect(markers.names, names(features))]
  }

  if (featuresType == "ENSEMBL_ID") {
    features <- subset(markersG, get("ENSEMBL ID") %in% rownames(expression))
    markers.names <- unique(features[, "Cell population"])
    features <- split(features[, "ENSEMBL ID"], features[, "Cell population"])
    missing.populations <- setdiff(markers.names, names(features))
    features <- features[intersect(markers.names, names(features))]
  }


  if (length(missing.populations) > 0) {
    warning(paste("Found no markers for population(s):", paste(missing.populations, collapse = ", ")))
  }
  t(appendSignatures(expression, features))
}

#' Test for Cell Population Infiltration
#'
#' Returns a matrix whose elements are p-values corresponding to the null hypothesis that samples are not infiltrated by the corresponding cell population.
#'
#' @param MCPcounterMatrix A matrix, usually a return from the MCPcounter.estimate method.
#' @param platform Expression platform used to produce the data. Supported are "133P2" (Affymetrix Human Genome 133 Plus 2.0), "133A" (Affymetrix Human Genome 133A), "HG1" (Affymetrix Human Gene 1.0ST). Other platforms are not supported. Data should ideally be log2-transformed and normalized with the fRMA bioconductor package. MCP-counter estimates from Affymetrix Human Genome 133 Plus 2.0 and 133A arrays should be computed using "affy_133P2_probesets" as identifiers, and "HUGO_symbols" or "ENTREZ_ID" for Affymetrix Human Gene 1.0ST.
#'
#' @return Matrix with samples in rows and cell populations in columns. Elements are p-values.
#' @author Etienne Becht
test_for_infiltration <- function(MCPcounterMatrix, platform = c("133P2", "133A", "HG1")[1]) {
  MCPcounterMatrix <- t(MCPcounterMatrix)
  params <- null_models[grep(platform, colnames(null_models))]
  rownames(params) <- null_models[, "Cell.population"]
  colnames(params) <- sub(platform, "", colnames(params), fixed = T)
  res <- sapply(colnames(MCPcounterMatrix), function(x) {
    pnorm(MCPcounterMatrix[, x], mean = params[x, "mu."], sd = params[x, "sigma."], lower.tail = F)
  })
  rownames(res) <- rownames(MCPcounterMatrix)
  res
}
