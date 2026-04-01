#' Identification of TME Cluster
#'
#' @description
#' Performs TME (Tumor Microenvironment) clustering analysis using various
#' clustering methods. Supports feature selection, scaling, and automatic
#' determination of optimal cluster number.
#'
#' @param input Data frame containing the input dataset.
#' @param features Vector of features to use for clustering. Default is NULL
#'   (uses all columns or pattern-selected columns).
#' @param pattern Regular expression pattern for selecting features.
#'   Default is NULL.
#' @param id Column name for identifiers. Default is NULL (uses row names).
#' @param method Clustering method. Default is "kmeans".
#' @param min_nc Minimum number of clusters to evaluate. Default is 2.
#' @param max.nc Maximum number of clusters to evaluate. Default is 6.
#' @param scale Logical indicating whether to scale features. Default is TRUE.
#'
#' @return Data frame with cluster assignments appended.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' tcga_stad_sig <- load_data("tcga_stad_sig")
#' res <- tme_cluster(
#'   input = tcga_stad_sig,
#'   pattern = "xCell",
#'   id = "ID",
#'   method = "kmeans"
#' )
#' }
tme_cluster <- function(input, features = NULL, pattern = NULL, id = NULL,
                        scale = TRUE, method = "kmeans", min_nc = 2, max.nc = 6) {
  input <- as.data.frame(input)

  # Extract IDs
  sample_ids <- if (is.null(id)) {
    rownames(input)
  } else {
    as.character(input[[id]])
  }

  # Select features
  if (is.null(features)) {
    if (!is.null(pattern)) {
      features <- colnames(input)[stringr::str_detect(colnames(input), pattern)]
    } else {
      features <- setdiff(colnames(input), id)
    }
  } else {
    features <- intersect(features, colnames(input))
  }

  if (length(features) == 0) {
    cli::cli_abort("No valid features selected for clustering")
  }

  input <- input[, features, drop = FALSE]

  # Scale features
  if (scale) {
    input <- as.data.frame(scale(input, scale = TRUE, center = TRUE))
    features <- feature_manipulation(data = input, feature = features, print_result = FALSE)
    input <- input[, features, drop = FALSE]
  }

  # Perform clustering
  rlang::check_installed("NbClust")
  res <- NbClust::NbClust(
    data = input,
    diss = NULL,
    distance = "euclidean",
    min.nc = min_nc,
    max.nc = max.nc,
    method = method,
    index = "kl",
    alphaBeale = 0.1
  )

  nc <- res$Best.nc
  cli::cli_alert_info("Best number of TME clusters: {nc[1]}")

  clusters <- res$Best.partition
  cli::cli_alert_info("Cluster distribution:")
  print(summary(as.factor(clusters)))

  out <- data.frame(ID = sample_ids, cluster = paste0("TME", clusters))
  input$ID <- sample_ids

  merge(out, input, by = "ID")
}
