#' The xCell Analysis Pipeline
#'
#' @description
#' Returns the xCell cell types enrichment scores for tumor microenvironment
#' deconvolution. Uses ssGSEA-based enrichment analysis with spillover
#' compensation to estimate cell type proportions from gene expression data.
#'
#' @param expr Gene expression data matrix with row names as gene symbols and
#'   columns as samples.
#' @param signatures GMT object of signatures. If `NULL`, uses xCell defaults.
#' @param genes Character vector of genes to use in the analysis. If `NULL`,
#'   uses xCell defaults.
#' @param spill Spillover object for adjusting scores. If `NULL`, uses xCell defaults.
#' @param rnaseq Logical indicating whether to use RNA-seq (TRUE) or array
#'   (FALSE) spillover parameters. Default is `TRUE`.
#' @param file.name Character string for saving scores. Default is `NULL`.
#' @param scale Logical indicating whether to use scaling. Default is `TRUE`.
#' @param alpha Numeric value to override spillover alpha parameter. Default is `0.5`.
#' @param save.raw Logical indicating whether to save raw scores. Default is `FALSE`.
#' @param cell.types.use Character vector of cell types to use. If `NULL`, uses
#'   all available cell types. Default is `NULL`.
#'
#' @return A matrix of adjusted xCell scores.
#'
#' @keywords internal
xCellAnalysis <- function(expr, signatures = NULL, genes = NULL,
                          spill = NULL, rnaseq = TRUE, file.name = NULL, scale = TRUE,
                          alpha = 0.5, save.raw = FALSE,
                          cell.types.use = NULL) {
  # Input validation
  if (is.null(expr)) {
    cli::cli_abort("{.arg expr} cannot be NULL")
  }
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    cli::cli_abort("{.arg expr} must be a matrix or data frame")
  }
  if (nrow(expr) == 0 || ncol(expr) == 0) {
    cli::cli_abort("{.arg expr} is empty")
  }
  if (is.null(rownames(expr))) {
    cli::cli_abort("{.arg expr} must have row names (gene symbols)")
  }
  if (!is.logical(rnaseq) || length(rnaseq) != 1) {
    cli::cli_abort("{.arg rnaseq} must be a single logical value")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    cli::cli_abort("{.arg alpha} must be a numeric value between 0 and 1")
  }

  rlang::check_installed("xCell")

  # Load default data if not provided
  signatures <- signatures %||% xCell::xCell.data$signatures
  genes <- genes %||% xCell::xCell.data$genes
  spill <- spill %||% if (rnaseq) {
    xCell::xCell.data$spill
  } else {
    xCell::xCell.data$spill.array
  }

  # Validate cell types if specified
  if (!is.null(cell.types.use)) {
    if (!is.character(cell.types.use)) {
      cli::cli_abort("{.arg cell.types.use} must be a character vector")
    }
    available_types <- rownames(spill$K)
    invalid_types <- setdiff(cell.types.use, available_types)
    if (length(invalid_types) > 0) {
      cli::cli_abort("Invalid cell types: {.val {head(invalid_types, 5)}}.
                     Available: {.val {head(available_types, 10)}}")
    }
  }

  fn <- if (is.null(file.name) || !save.raw) {
    NULL
  } else {
    paste0(file.name, "_RAW.txt")
  }

  scores <- rawEnrichmentAnalysis(expr, signatures, genes, fn)
  scores.transformed <- transformScores(scores, spill$fv, scale)

  if (is.null(cell.types.use)) {
    scores.adjusted <- spillOver(scores.transformed, spill$K, alpha, file.name)
    scores.adjusted <- microenvironmentScores(scores.adjusted)
  } else {
    scores.adjusted <- spillOver(scores.transformed[cell.types.use, ],
                                  spill$K, alpha, file.name)
  }

  scores.adjusted
}

#' Calculate Raw xCell Enrichment Scores
#'
#' @description
#' Returns the raw xCell cell types enrichment scores using ssGSEA.
#'
#' @param expr Gene expression data matrix with row names as gene symbols.
#' @param signatures GMT object of signatures.
#' @param genes Character vector of genes to use.
#' @param file.name Character string for saving scores. Default is `NULL`.
#'
#' @return Matrix of raw xCell scores.
#'
#' @keywords internal
rawEnrichmentAnalysis <- function(expr, signatures, genes, file.name = NULL) {
  # Reduce expression dataset to required genes
  shared.genes <- intersect(rownames(expr), genes)
  cli::cli_alert_info("Number of genes: {length(shared.genes)}")

  if (length(shared.genes) < 5000) {
    cli::cli_abort("Not enough genes. Need at least 5000, found {length(shared.genes)}")
  }

  expr <- expr[shared.genes, ]
  expr <- apply(expr, 2, rank)

  # Run ssGSEA analysis
  gsva_info <- gsva_use_new_api()
  use_new_api <- gsva_info$use_new_api

  if (use_new_api) {
    params <- GSVA::gsvaParam(
      as.matrix(expr),
      signatures,
      minSize = 1,
      maxSize = Inf,
      kcdf = "Gaussian",
      tau = 1,
      maxDiff = TRUE,
      absRanking = FALSE
    )

    rlang::check_installed("BiocParallel")
    scores <- GSVA::gsva(
      params,
      verbose = TRUE,
      BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
    )
  } else {
    scores <- GSVA::gsva(
      expr,
      signatures,
      method = "ssgsea",
      ssgsea.norm = FALSE
    )
  }

  scores <- scores - apply(scores, 1, min)

  # Combine signatures for same cell types
  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]
  agg <- stats::aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]

  if (!is.null(file.name)) {
    utils::write.table(scores,
      file = file.name, sep = "\t",
      col.names = NA, quote = FALSE
    )
  }

  scores
}

#' Transform Raw Scores to Fractions
#'
#' @description
#' Transforms raw xCell scores to estimated cell fractions.
#'
#' @param scores Raw scores from rawEnrichmentAnalysis.
#' @param fit.vals Calibration values from spill object.
#' @param scale Logical indicating whether to use scaling. Default is `TRUE`.
#' @param fn Character string for saving scores. Default is `NULL`.
#'
#' @return Transformed xCell scores matrix.
#'
#' @keywords internal
transformScores <- function(scores, fit.vals, scale = TRUE, fn = NULL) {
  rows <- rownames(scores)[rownames(scores) %in% rownames(fit.vals)]
  tscores <- scores[rows, ]
  minX <- apply(tscores, 1, min)
  A <- rownames(tscores)
  tscores <- (as.matrix(tscores) - minX) / 5000
  tscores[tscores < 0] <- 0

  if (!scale) {
    fit.vals[A, 3] <- 1
  }

  tscores <- (tscores^fit.vals[A, 2]) / (fit.vals[A, 3] * 2)

  if (!is.null(fn)) {
    utils::write.table(format(tscores, digits = 4),
      file = fn, sep = "\t",
      col.names = NA, quote = FALSE
    )
  }

  tscores
}

#' Adjust Scores Using Spillover Compensation
#'
#' @description
#' Adjusts xCell scores using spillover compensation matrix.
#'
#' @param transformedScores Transformed scores from transformScores.
#' @param K Spillover matrix.
#' @param alpha Spillover alpha parameter. Default is `0.5`.
#' @param file.name Character string for saving scores. Default is `NULL`.
#'
#' @return Adjusted xCell scores matrix.
#'
#' @keywords internal
spillOver <- function(transformedScores, K, alpha = 0.5, file.name = NULL) {
  K <- K * alpha
  diag(K) <- 1
  rows <- rownames(transformedScores)[rownames(transformedScores) %in% rownames(K)]

  rlang::check_installed("limSolve")

  scores <- apply(transformedScores[rows, , drop = FALSE], 2, function(x) {
    G <- diag(nrow(K[rows, rows]))
    H <- rep(0, nrow(G))
    res <- limSolve::lsei(
      A = K[rows, rows],
      B = x,
      G = G,
      H = H,
      verbose = FALSE
    )
    pmax(res$X, 0)
  })

  scores[scores < 0] <- 0
  rownames(scores) <- rows

  if (!is.null(file.name)) {
    scores <- round(scores * 10000) / 10000
    utils::write.table(scores,
      file = file.name, sep = "\t",
      col.names = NA, quote = FALSE
    )
  }

  scores
}

#' Calculate Microenvironment Scores
#'
#' @description
#' Calculates combined immune and stroma scores from adjusted xCell scores.
#'
#' @param adjustedScores Adjusted xCell scores matrix.
#'
#' @return Matrix with additional microenvironment score rows.
#'
#' @keywords internal
microenvironmentScores <- function(adjustedScores) {
  immune_cells <- c(
    "B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC",
    "Eosinophils", "Macrophages", "Monocytes", "Mast cells",
    "Neutrophils", "NK cells"
  )
  stroma_cells <- c("Adipocytes", "Endothelial cells", "Fibroblasts")

  ImmuneScore <- apply(adjustedScores[immune_cells, ], 2, sum) / 1.5
  StromaScore <- apply(adjustedScores[stroma_cells, ], 2, sum) / 2
  MicroenvironmentScore <- ImmuneScore + StromaScore

  rbind(adjustedScores,
    ImmuneScore = ImmuneScore,
    StromaScore = StromaScore, MicroenvironmentScore = MicroenvironmentScore
  )
}

#' Calculate Significance P-values Using Beta Distribution
#'
#' @description
#' Calculates FDR-adjusted p-values for the null hypothesis that a cell type
#' is not present in the mixture.
#'
#' @param scores xCell scores matrix.
#' @param beta_params Pre-calculated beta distribution parameters.
#' @param rnaseq Logical for RNA-seq vs array parameters. Default is `TRUE`.
#' @param file.name Character string for saving p-values. Default is `NULL`.
#'
#' @return Matrix of p-values.
#'
#' @keywords internal
xCellSignifcanceBetaDist <- function(scores, beta_params = NULL, rnaseq = TRUE,
                                     file.name = NULL) {
  rlang::check_installed("xCell")

  beta_params <- beta_params %||% if (rnaseq) {
    xCell::xCell.data$spill$beta_params
  } else {
    xCell::xCell.data$spill.array$beta_params
  }

  scores <- scores[rownames(scores) %in%
    colnames(xCell::xCell.data$spill$beta_params[[1]]), ]
  pvals <- matrix(0, nrow(scores), ncol(scores))
  rownames(pvals) <- rownames(scores)
  eps <- 1e-3

  for (i in seq_len(nrow(scores))) {
    ct <- rownames(scores)[i]
    beta_dist <- c()

    for (bp in beta_params) {
      if (sum(bp[, i] == 0)) {
        bd <- matrix(eps, 1, 100000)
      } else {
        bd <- stats::rbeta(100000, bp[1, ct], bp[2, ct])
        bd <- ((1 + eps) * (bp[3, ct])) * bd
      }
      beta_dist <- c(beta_dist, bd)
    }

    pvals[i, ] <- 1 - mapply(scores[i, ], FUN = function(x) mean(x > beta_dist))
  }

  if (!is.null(file.name)) {
    utils::write.table(pvals,
      file = file.name, quote = FALSE,
      row.names = TRUE, sep = "\t", col.names = NA
    )
  }

  pvals
}

#' Calculate Significance Using Random Matrix
#'
#' @description
#' Calculates FDR-adjusted p-values using a random shuffled matrix.
#'
#' @param scores xCell scores matrix.
#' @param expr Input expression matrix.
#' @param spill Spillover object.
#' @param alpha Spillover alpha parameter. Default is `0.5`.
#' @param nperm Number of permutations. Default is `250`.
#' @param file.name Character string for saving p-values. Default is `NULL`.
#'
#' @return List containing p-values, shuffled xCell scores, shuffled expression,
#'   and beta distributions.
#'
#' @keywords internal
xCellSignifcanceRandomMatrix <- function(scores, expr, spill, alpha = 0.5,
                                         nperm = 250, file.name = NULL) {
  rlang::check_installed("xCell")

  shuff_expr <- mapply(seq_len(nperm),
                       FUN = function(x) sample(nrow(expr), nrow(expr)))
  rownames(shuff_expr) <- sample(rownames(expr))
  shuff_xcell <- xCellAnalysis(shuff_expr, spill = spill, alpha = alpha)
  shuff_xcell <- shuff_xcell[rownames(scores), ]

  pvals <- matrix(0, nrow(scores), ncol(scores))
  beta_dist <- matrix(0, nrow(scores), 100000)
  eps <- 1e-3

  for (i in seq_len(nrow(scores))) {
    x <- shuff_xcell[i, ]
    if (stats::sd(x) < eps) {
      beta_dist[i, ] <- rep(eps, 100000)
    } else {
      x <- x + eps
      rlang::check_installed("MASS")
      beta_params <- MASS::fitdistr(x / ((1 + 2 * eps) * (max(x))) + eps,
        "beta", list(shape1 = 1, shape2 = 1),
        lower = eps
      )
      beta_dist[i, ] <- stats::rbeta(
        100000, beta_params$estimate[1],
        beta_params$estimate[2]
      )
      beta_dist[i, ] <- ((1 + 2 * eps) * (max(x))) * beta_dist[i, ]
    }

    pvals[i, ] <- 1 - unlist(lapply(scores[i, ],
      FUN = function(x) mean(x > beta_dist[i, ])
    ))
  }

  rownames(pvals) <- rownames(scores)
  colnames(pvals) <- colnames(scores)
  rownames(shuff_xcell) <- rownames(scores)
  rownames(beta_dist) <- rownames(scores)

  if (!is.null(file.name)) {
    utils::write.table(pvals,
      file = file.name, quote = FALSE,
      row.names = TRUE, sep = "\t", col.names = NA
    )
  }

  list(
    pvals = pvals, shuff_xcell = shuff_xcell,
    shuff_expr = shuff_expr, beta_dist = beta_dist
  )
}
