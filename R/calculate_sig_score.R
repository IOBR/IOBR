#' Signature Score Calculation Methods
#'
#' @description
#' A named vector of supported methods for calculating signature scores.
#'
#' @format Named character vector:
#' \describe{
#'   \item{PCA}{Principal Component Analysis method ("pca")}
#'   \item{ssGSEA}{Single-sample Gene Set Enrichment Analysis ("ssgsea")}
#'   \item{z-score}{Z-score transformation method ("zscore")}
#'   \item{Integration}{Integration of multiple methods ("integration")}
#' }
#'
#' @export
#'
#' @examples
#' signature_score_calculation_methods
#' signature_score_calculation_methods["PCA"]
signature_score_calculation_methods <- c(
  "PCA"         = "pca",
  "ssGSEA"      = "ssgsea",
  "z-score"     = "zscore",
  "Integration" = "integration"
)

# Helper: Prepare pdata and match with eset
.prepare_pdata <- function(pdata, eset, column_of_sample) {
  if (is.null(pdata)) {
    pdata <- data.frame(
      Index = seq_len(ncol(eset)),
      ID = colnames(eset),
      stringsAsFactors = FALSE
    )
  } else {
    pdata <- as.data.frame(pdata)

    # Rename existing ID column to avoid conflict
    if ("ID" %in% colnames(pdata) && column_of_sample != "ID") {
      colnames(pdata)[colnames(pdata) == "ID"] <- "ID2"
      cli::cli_alert_warning(
        "Renamed existing 'ID' column to 'ID2' to avoid conflicts"
      )
    }

    if (column_of_sample %in% colnames(pdata)) {
      colnames(pdata)[colnames(pdata) == column_of_sample] <- "ID"
    }
  }

  # Match pdata with eset
  pdata <- pdata[pdata$ID %in% colnames(eset), , drop = FALSE]
  eset <- eset[, colnames(eset) %in% pdata$ID, drop = FALSE]
  eset <- eset[, match(pdata$ID, colnames(eset)), drop = FALSE]

  list(pdata = pdata, eset = eset)
}

# Helper: Filter signatures by gene count
.filter_signatures <- function(signature, eset, mini_gene_count) {
  gene_counts <- vapply(
    signature,
    function(x) sum(x %in% rownames(eset)),
    integer(1)
  )
  signature[gene_counts >= mini_gene_count]
}

# Helper: Add TME combined scores
.add_tme_scores <- function(pdata, suffix = "") {
  pca_suffix <- ifelse(nzchar(suffix), paste0("_", suffix), "")

  if (all(c(
    paste0("TMEscoreA_CIR", pca_suffix),
    paste0("TMEscoreB_CIR", pca_suffix)
  ) %in% colnames(pdata))) {
    pdata[[paste0("TMEscore_CIR", pca_suffix)]] <-
      pdata[[paste0("TMEscoreA_CIR", pca_suffix)]] -
      pdata[[paste0("TMEscoreB_CIR", pca_suffix)]]
  }

  if (all(c(
    paste0("TMEscoreA_plus", pca_suffix),
    paste0("TMEscoreB_plus", pca_suffix)
  ) %in% colnames(pdata))) {
    pdata[[paste0("TMEscore_plus", pca_suffix)]] <-
      pdata[[paste0("TMEscoreA_plus", pca_suffix)]] -
      pdata[[paste0("TMEscoreB_plus", pca_suffix)]]
  }

  pdata
}

#' Calculate Signature Score Using PCA Method
#'
#' @description
#' Computes signature scores using Principal Component Analysis.
#' The first principal component is used as the signature score.
#'
#' @param pdata Data frame with phenotype data. If `NULL`, created from
#'   `eset` column names.
#' @param eset Expression matrix (genes as rows, samples as columns).
#' @param signature List of gene signatures.
#' @param mini_gene_count Minimum genes required per signature. Default is 3.
#' @param column_of_sample Column in `pdata` with sample IDs. Default is `"ID"`.
#' @param adjust_eset Logical: remove problematic features. Default is `FALSE`.
#'
#' @return Tibble with signature scores.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' set.seed(123)
#' eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(eset) <- paste0("Gene", 1:100)
#' colnames(eset) <- paste0("Sample", 1:10)
#' signature <- list(
#'   Signature1 = paste0("Gene", 1:10),
#'   Signature2 = paste0("Gene", 11:20)
#' )
#' \donttest{
#' result <- calculate_sig_score_pca(eset = eset, signature = signature)
#' }
calculate_sig_score_pca <- function(pdata = NULL,
                                    eset,
                                    signature,
                                    mini_gene_count = 3,
                                    column_of_sample = "ID",
                                    adjust_eset = FALSE) {
  cli::cli_alert_info("Calculating signature scores using PCA method")

  # Prepare data
  prepared <- .prepare_pdata(pdata, eset, column_of_sample)
  pdata <- prepared$pdata
  eset <- prepared$eset

  if (ncol(eset) == 0) {
    cli::cli_abort("No matching samples between pdata and eset")
  }

  # Preprocessing for small datasets
  if (ncol(eset) < 5000) {
    eset <- log2eset(eset)
    check_eset(eset)
  }

  if (adjust_eset) {
    feas <- feature_manipulation(data = eset, is_matrix = TRUE)
    eset <- eset[rownames(eset) %in% feas, , drop = FALSE]
  }

  # Filter signatures
  signature <- .filter_signatures(signature, eset, max(mini_gene_count, 2))

  if (length(signature) == 0) {
    cli::cli_warn("No signatures have enough genes. Returning pdata unchanged.")
    return(tibble::as_tibble(pdata))
  }

  cli::cli_alert_info("Calculating scores for {length(signature)} signature(s)")

  # Calculate scores
  sig_scores <- vapply(names(signature), function(sig) {
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    if (length(genes) == 0) {
      return(rep(NA_real_, ncol(eset)))
    }
    sigScore(eset[genes, , drop = FALSE], methods = "PCA")
  }, numeric(ncol(eset)))

  # Add to pdata
  pdata <- cbind(pdata, as.data.frame(sig_scores))
  pdata <- .add_tme_scores(pdata)

  if ("Index" %in% colnames(pdata)) {
    pdata <- pdata[, colnames(pdata) != "Index", drop = FALSE]
  }

  tibble::as_tibble(pdata)
}

#' Calculate Signature Score Using Z-Score Method
#'
#' @description
#' Computes signature scores using z-score transformation.
#'
#' @inheritParams calculate_sig_score_pca
#'
#' @return Tibble with signature scores.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' set.seed(123)
#' eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(eset) <- paste0("Gene", 1:100)
#' colnames(eset) <- paste0("Sample", 1:10)
#' signature <- list(
#'   Signature1 = paste0("Gene", 1:10),
#'   Signature2 = paste0("Gene", 11:20)
#' )
#' \donttest{
#' result <- calculate_sig_score_zscore(eset = eset, signature = signature)
#' }
calculate_sig_score_zscore <- function(pdata = NULL,
                                       eset,
                                       signature,
                                       mini_gene_count = 3,
                                       column_of_sample = "ID",
                                       adjust_eset = FALSE) {
  cli::cli_alert_info("Calculating signature scores using z-score method")

  # Prepare data
  prepared <- .prepare_pdata(pdata, eset, column_of_sample)
  pdata <- prepared$pdata
  eset <- prepared$eset

  if (ncol(eset) < 5000) {
    eset <- log2eset(eset)
    check_eset(eset)
  }

  if (adjust_eset) {
    feas <- feature_manipulation(data = eset, is_matrix = TRUE)
    eset <- eset[rownames(eset) %in% feas, , drop = FALSE]
  }

  # Filter signatures
  signature <- .filter_signatures(signature, eset, max(mini_gene_count, 2))

  if (length(signature) == 0) {
    cli::cli_warn("No signatures have enough genes. Returning pdata unchanged.")
    return(tibble::as_tibble(pdata))
  }

  cli::cli_alert_info("Calculating scores for {length(signature)} signature(s)")

  # Calculate scores
  sig_scores <- vapply(names(signature), function(sig) {
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    if (length(genes) == 0) {
      return(rep(NA_real_, ncol(eset)))
    }
    sigScore(eset[genes, , drop = FALSE], methods = "zscore")
  }, numeric(ncol(eset)))

  pdata <- cbind(pdata, as.data.frame(sig_scores))
  pdata <- .add_tme_scores(pdata)

  if ("Index" %in% colnames(pdata)) {
    pdata <- pdata[, colnames(pdata) != "Index", drop = FALSE]
  }

  tibble::as_tibble(pdata)
}

#' Calculate Signature Score Using ssGSEA Method
#'
#' @description
#' Computes signature scores using single-sample Gene Set Enrichment Analysis.
#'
#' @inheritParams calculate_sig_score_pca
#' @param parallel.size Number of parallel workers. Default is 1.
#'
#' @return Tibble with signature scores.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' set.seed(123)
#' eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(eset) <- paste0("Gene", 1:100)
#' colnames(eset) <- paste0("Sample", 1:10)
#' signature <- list(
#'   Signature1 = paste0("Gene", 1:15),
#'   Signature2 = paste0("Gene", 16:30)
#' )
#' \donttest{
#' result <- calculate_sig_score_ssgsea(eset = eset, signature = signature)
#' }
calculate_sig_score_ssgsea <- function(pdata = NULL,
                                       eset,
                                       signature,
                                       mini_gene_count = 5,
                                       column_of_sample = "ID",
                                       adjust_eset = FALSE,
                                       parallel.size = 1L) {
  cli::cli_alert_info("Calculating signature scores using ssGSEA method")

  # ssGSEA needs more genes for robustness
  min_genes <- max(mini_gene_count, 5)

  # Filter signatures early
  signature <- .filter_signatures(signature, eset, min_genes)

  if (length(signature) == 0) {
    cli::cli_warn(paste0(
      "No signatures have enough genes (min: ",
      min_genes, "). Returning pdata unchanged."
    ))
    return(tibble::as_tibble(pdata))
  }

  # Prepare data
  prepared <- .prepare_pdata(pdata, eset, column_of_sample)
  pdata <- prepared$pdata
  eset <- prepared$eset

  if (ncol(eset) < 5000) {
    eset <- log2eset(eset)
    check_eset(eset)
  }

  if (adjust_eset) {
    feas <- feature_manipulation(data = eset, is_matrix = TRUE)
    eset <- eset[rownames(eset) %in% feas, , drop = FALSE]
  }

  cli::cli_alert_info("Calculating scores for {length(signature)} signature(s)")

  # Run ssGSEA with appropriate API
  rlang::check_installed("GSVA")
  gsva_info <- gsva_use_new_api()

  if (gsva_info$use_new_api) {
    rlang::check_installed("BiocParallel")

    bp <- if (parallel.size > 1) {
      if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(workers = parallel.size, progressbar = TRUE)
      } else {
        BiocParallel::MulticoreParam(
          workers = parallel.size,
          progressbar = TRUE
        )
      }
    } else {
      BiocParallel::SerialParam(progressbar = TRUE)
    }

    param <- GSVA::ssgseaParam(
      exprData = as.matrix(eset),
      geneSets = signature,
      minSize = min_genes,
      maxSize = Inf,
      normalize = TRUE
    )
    res <- GSVA::gsva(param, verbose = TRUE, BPPARAM = bp)
  } else {
    res <- GSVA::gsva(
      as.matrix(eset),
      signature,
      method = "ssgsea",
      kcdf = "Gaussian",
      min.sz = min_genes,
      ssgsea.norm = TRUE,
      parallel.sz = parallel.size
    )
  }

  # Process results
  res <- as.data.frame(t(res))
  res <- .add_tme_scores(res)
  res <- tibble::rownames_to_column(res, var = "ID")

  pdata <- merge(pdata, res, by = "ID", all = FALSE)

  if ("Index" %in% colnames(pdata)) {
    pdata <- pdata[, colnames(pdata) != "Index", drop = FALSE]
  }

  tibble::as_tibble(pdata)
}

#' Calculate Signature Score Using Integration Method
#'
#' @description
#' Computes signature scores using PCA, z-score, and ssGSEA methods combined.
#'
#' @inheritParams calculate_sig_score_ssgsea
#'
#' @return Tibble with signature scores from all three methods.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' set.seed(123)
#' eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(eset) <- paste0("Gene", 1:100)
#' colnames(eset) <- paste0("Sample", 1:10)
#' signature <- list(
#'   Signature1 = paste0("Gene", 1:15),
#'   Signature2 = paste0("Gene", 16:30)
#' )
#' \donttest{
#' result <- calculate_sig_score_integration(eset = eset, signature = signature)
#' }
calculate_sig_score_integration <- function(pdata = NULL,
                                            eset,
                                            signature,
                                            mini_gene_count = 2,
                                            column_of_sample = "ID",
                                            adjust_eset = FALSE,
                                            parallel.size = 1L) {
  cli::cli_alert_info("Calculating signature scores using PCA, z-score, and ssGSEA methods")

  # Filter signatures
  signature <- .filter_signatures(signature, eset, mini_gene_count)

  if (length(signature) == 0) {
    cli::cli_abort("No signatures have enough genes (min: {mini_gene_count})")
  }

  # Prepare data
  prepared <- .prepare_pdata(pdata, eset, column_of_sample)
  pdata <- prepared$pdata
  eset <- prepared$eset

  if (ncol(eset) < 5000) {
    eset <- log2eset(eset)
    check_eset(eset)
  }

  if (adjust_eset) {
    feas <- feature_manipulation(data = eset, is_matrix = TRUE)
    eset <- eset[rownames(eset) %in% feas, , drop = FALSE]
  }

  goi <- names(signature)

  # Step 1: PCA
  cli::cli_alert_info("Step 1/3: PCA method")
  for (sig in goi) {
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    pdata[[paste0(sig, "_PCA")]] <- sigScore(eset[genes, , drop = FALSE], methods = "PCA")
  }
  pdata <- .add_tme_scores(pdata, "PCA")

  # Step 2: z-score
  cli::cli_alert_info("Step 2/3: z-score method")
  for (sig in goi) {
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    pdata[[paste0(sig, "_zscore")]] <- sigScore(eset[genes, , drop = FALSE],
      methods = "zscore"
    )
  }
  pdata <- .add_tme_scores(pdata, "zscore")

  # Step 3: ssGSEA
  cli::cli_alert_info("Step 3/3: ssGSEA method")
  ssgsea_min <- max(mini_gene_count, 5)
  sig_ssgsea <- signature[vapply(
    signature,
    function(x) sum(x %in% rownames(eset)),
    integer(1)
  ) >= ssgsea_min]

  if (length(sig_ssgsea) > 0) {
    rlang::check_installed("GSVA")
    gsva_info <- gsva_use_new_api()

    if (gsva_info$use_new_api) {
      rlang::check_installed("BiocParallel")
      param <- GSVA::ssgseaParam(
        exprData = as.matrix(eset),
        geneSets = sig_ssgsea,
        minSize = ssgsea_min,
        maxSize = Inf,
        normalize = TRUE
      )
      res <- GSVA::gsva(param,
        verbose = TRUE,
        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
      )
    } else {
      res <- GSVA::gsva(
        as.matrix(eset),
        sig_ssgsea,
        method = "ssgsea",
        kcdf = "Gaussian",
        min.sz = ssgsea_min,
        ssgsea.norm = TRUE
      )
    }

    res <- as.data.frame(t(res))
    res <- .add_tme_scores(res)
    colnames(res) <- paste0(colnames(res), "_ssGSEA")
    res <- tibble::rownames_to_column(res, var = "ID")
    pdata <- merge(pdata, res, by = "ID", all = FALSE)
  }

  tibble::as_tibble(pdata)
}

#' Calculate Signature Score
#'
#' @description
#' Main interface for calculating signature scores from gene expression data.
#' Supports PCA, z-score, ssGSEA, and integration methods.
#'
#' @param pdata Phenotype data. If `NULL`, created from `eset` column names.
#' @param eset Expression matrix (CPM, TPM, RPKM, FPKM, etc.).
#' @param signature List of gene signatures. Can also be a character string
#'   naming a built-in signature collection (e.g., `"signature_collection"`,
#'   `"signature_tme"`, `"go_bp"`, `"kegg"`, `"hallmark"`).
#' @param method Scoring method: `"pca"`, `"ssgsea"`, `"zscore"`,
#'   or `"integration"`. Default is `"pca"`.
#' @param mini_gene_count Minimum genes required per signature. Default is 3
#'   (or 5 for ssGSEA).
#' @param column_of_sample Column with sample IDs in `pdata`. Default is `"ID"`.
#' @param print_gene_proportion Logical: print gene coverage. Default is
#'   `FALSE`.
#' @param print_filtered_signatures Logical: print filtered signatures.
#'   Default is `FALSE`.
#' @param adjust_eset Logical: clean problematic features. Default is `FALSE`.
#' @param parallel.size Parallel workers for ssGSEA. Default is 1.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return Tibble with phenotype data and signature scores.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @references
#' \enumerate{
#'   \item Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation
#'     analysis. BMC Bioinformatics. 2013;14:7.
#'   \item Mariathasan S, et al. TGF\eqn{\beta} attenuates tumour response to PD-L1 blockade.
#'     Nature. 2018;554:544-548.
#' }
#'
#' @examples
#' set.seed(123)
#' eset <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(eset) <- paste0("Gene", 1:100)
#' colnames(eset) <- paste0("Sample", 1:10)
#' signature <- list(
#'   Signature1 = paste0("Gene", 1:10),
#'   Signature2 = paste0("Gene", 11:20)
#' )
#' \donttest{
#' result <- calculate_sig_score(eset = eset, signature = signature, method = "pca")
#' }
calculate_sig_score <- function(pdata = NULL,
                                eset,
                                signature = NULL,
                                method = c(
                                  "pca", "ssgsea", "zscore",
                                  "integration"
                                ),
                                mini_gene_count = 3,
                                column_of_sample = "ID",
                                print_gene_proportion = FALSE,
                                print_filtered_signatures = FALSE,
                                adjust_eset = FALSE,
                                parallel.size = 1L,
                                ...) {
  # Validate signature
  if (is.null(signature)) {
    cli::cli_abort(c(
      "Please provide a signature list.",
      "i" = "Available: signature_tme, signature_collection, go_bp, \
              kegg, hallmark"
    ))
  }

  # Load signature if character
  if (is.character(signature) && length(signature) == 1) {
    builtin_sigs <- c(
      "signature_collection", "signature_tme", "go_bp",
      "kegg", "hallmark"
    )
    if (signature %in% builtin_sigs) {
      signature <- load_data(signature)
    }
  }

  # Clean signature list
  if (is.list(signature)) {
    signature <- lapply(signature, function(x) {
      x <- stats::na.omit(as.character(x))
      unique(x[nzchar(x)])
    })
  }

  # Print gene proportion if requested
  if (print_gene_proportion) {
    coverage <- lapply(signature, function(x) {
      prop <- mean(x %in% rownames(eset))
      cli::cli_alert_info("{length(x)} genes, {prop*100:.1f}% in expression data")
    })
  }

  # Print filtered signatures if requested
  if (print_filtered_signatures) {
    gene_counts <- vapply(signature, function(x) sum(x %in% rownames(eset)), integer(1))
    filtered <- signature[gene_counts <= mini_gene_count]

    cli::cli_alert_info("{length(filtered)} signature(s) filtered (<{mini_gene_count} genes)")
    if (length(filtered) > 0 && length(filtered) <= 10) {
      cli::cli_text("Filtered: {.val {names(filtered)}}")
    } else if (length(filtered) > 10) {
      cli::cli_text("First 10 filtered: {.val {names(filtered)[1:10]}}")
    }
  }

  # Validate and dispatch
  method <- rlang::arg_match(method)

  switch(method,
    pca = calculate_sig_score_pca(
      pdata, eset,
      signature = signature,
      mini_gene_count = mini_gene_count,
      column_of_sample = column_of_sample,
      adjust_eset = adjust_eset, ...
    ),
    ssgsea = calculate_sig_score_ssgsea(
      pdata, eset,
      signature = signature,
      mini_gene_count = mini_gene_count,
      column_of_sample = column_of_sample,
      adjust_eset = adjust_eset,
      parallel.size = parallel.size, ...
    ),
    zscore = calculate_sig_score_zscore(
      pdata, eset,
      signature = signature,
      mini_gene_count = mini_gene_count,
      column_of_sample = column_of_sample,
      adjust_eset = adjust_eset, ...
    ),
    integration = calculate_sig_score_integration(
      pdata, eset,
      signature = signature,
      mini_gene_count = mini_gene_count,
      column_of_sample = column_of_sample,
      adjust_eset = adjust_eset,
      parallel.size = parallel.size, ...
    )
  )
}
