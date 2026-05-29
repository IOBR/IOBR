#' TME Deconvolution Methods
#'
#' @description
#' A named vector of supported tumor microenvironment (TME) deconvolution
#' methods in the IOBR package.
#'
#' @details
#' The methods currently supported are:
#' \itemize{
#'   \item `mcpcounter`: MCP-counter for immune and stromal cell populations
#'   \item `epic`: EPIC for immune, stromal, and cancer cell fractions
#'   \item `xcell`: xCell for 64 immune and stromal cell types
#'   \item `cibersort`: CIBERSORT for 22 immune cell types
#'   \item `cibersort_abs`: CIBERSORT in absolute mode
#'   \item `ips`: Immunophenoscore calculation
#'   \item `estimate`: ESTIMATE for stromal/immune/estimate scores
#'   \item `svr`: Support Vector Regression (custom reference)
#'   \item `lsei`: Least Squares with Equality/Inequality constraints
#'   \item `timer`: TIMER for cancer-specific immune estimation
#'   \item `quantiseq`: quanTIseq for RNA-seq immune deconvolution
#' }
#'
#' @format A named character vector where names are display names and values
#'   are internal method names.
#'
#' @return Named character vector of available deconvolution methods.
#'
#' @export
#'
#' @examples
#' # List all available TME deconvolution methods
#' tme_deconvolution_methods
#'
#' # Get internal method name for a specific method
#' tme_deconvolution_methods["MCPcounter"]
tme_deconvolution_methods <- c(
  "MCPcounter"         = "mcpcounter",
  "EPIC"               = "epic",
  "xCell"              = "xcell",
  "CIBERSORT"          = "cibersort",
  "CIBERSORT Absolute" = "cibersort_abs",
  "IPS"                = "ips",
  "ESTIMATE"           = "estimate",
  "SVR"                = "svr",
  "lsei"               = "lsei",
  "TIMER"              = "timer",
  "quanTIseq"          = "quantiseq"
)

# Helper: Add project ID and format result
.format_deconv_result <- function(res, project, suffix) {
  if (is.null(res) || ncol(res) == 0 || nrow(res) == 0) return(NULL)

  colnames(res) <- gsub("\\.", "_", colnames(res))
  colnames(res) <- gsub(" ", "_", colnames(res))
  colnames(res) <- paste0(colnames(res), "_", suffix)

  if (!is.null(project)) {
    res$ProjectID <- project
    res <- res[, c(ncol(res), seq_len(ncol(res) - 1)), drop = FALSE]
  }

  tibble::rownames_to_column(as.data.frame(res), var = "ID")
}

#' Deconvolve Immune Microenvironment Using xCell
#'
#' @description
#' Estimates immune cell fractions using the xCell algorithm. xCell provides
#' cell type enrichment scores for 64 immune and stromal cell types from gene
#' expression data.
#'
#' @param eset Gene expression matrix with HGNC gene symbols as row names and
#'   samples as columns.
#' @param project Optional project name to add as `ProjectID` column.
#'   Default is `NULL`.
#' @param arrays Logical indicating microarray data (`TRUE`) or RNA-seq
#'   (`FALSE`). Default is `FALSE`.
#'
#' @return Data frame with xCell enrichment scores. Cell type columns are
#'   suffixed with `_xCell`.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' xcell <- load_data("xCell.data")
#' if (!is.null(xcell)) {
#'   set.seed(123)
#'   sim_eset <- matrix(rnorm(length(xcell$genes) * 3), length(xcell$genes), 3)
#'   rownames(sim_eset) <- xcell$genes
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'   
#'   # Run deconvolution
#'   result <- deconvo_xcell(eset = sim_eset, project = "TCGA-STAD")
#'   if (!is.null(result)) head(result)
#' }
deconvo_xcell <- function(eset, project = NULL, arrays = FALSE) {
  if (is.null(eset)) return(NULL)
  cli::cli_alert_info("Running xCell deconvolution")

  # Validate row names
  rn <- rownames(eset)
  if (is.null(rn) || length(rn) == 0) {
    cli::cli_abort(c(
      "xCell requires gene symbols as row names.",
      "i" = "Your input has no row names.",
      "*" = "For RNA-seq: use HGNC gene symbols as row names with arrays=FALSE",
      "*" = "For microarray: set arrays=TRUE or map probes to gene symbols"
    ))
  }

  # Check identifier format
  ensg_frac <- mean(grepl("^ENSG", rn))
  ensver_frac <- mean(grepl("\\.\\d+$", rn))
  probe_frac <- mean(grepl("(_at$|_s_at$|_x_at$|^AFFX)", rn,
    ignore.case = TRUE
  ))

  if (ensg_frac > 0.5 || probe_frac > 0.5) {
    cli::cli_abort(c(
      "Gene identifier format appears incompatible with xCell.",
      "i" = "Detected: {ensg_frac*100:.1f}% Ensembl-like, \
              {probe_frac*100:.1f}% probe-like",
      "*" = "xCell requires HGNC gene symbols for RNA-seq mode",
      "*" = "Use anno_eset() to convert identifiers to gene symbols"
    ))
  }

  if (any(grepl("\\s", rn))) {
    cli::cli_warn("Row names contain whitespace; consider trimws()")
  }
  if (all(rn == tolower(rn))) {
    cli::cli_warn("Row names are lowercase; xCell uses uppercase HGNC symbols")
  }

  # Run xCell
  rnaseq <- !arrays
  res <- xCellAnalysis(eset, rnaseq = rnaseq)
  if (is.null(res)) return(NULL)
  res <- as.data.frame(t(res))

  .format_deconv_result(res, project, "xCell")
}

#' Deconvolve Immune Microenvironment Using MCP-Counter
#'
#' @description
#' Estimates immune cell abundances using MCP-counter.
#'
#' @param eset Gene expression matrix with HGNC symbols as row names.
#' @param project Optional project name. Default is `NULL`.
#'
#' @return Data frame with MCP-counter scores. Columns suffixed with
#'   `_MCPcounter`.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' mcp_genes <- load_data("mcp_genes")
#' if (!is.null(mcp_genes)) {
#'   set.seed(123)
#'   sim_eset <- matrix(rnorm(nrow(mcp_genes) * 3), nrow(mcp_genes), 3)
#'   rownames(sim_eset) <- mcp_genes$`HUGO symbols`
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'   
#'   # Run deconvolution
#'   result <- deconvo_mcpcounter(eset = sim_eset, project = "TCGA-STAD")
#'   if (!is.null(result)) head(result)
#' }
deconvo_mcpcounter <- function(eset, project = NULL) {
  if (is.null(eset)) return(NULL)
  cli::cli_alert_info("Running MCP-counter deconvolution")

  mcp_genes <- load_data("mcp_genes")
  if (is.null(mcp_genes)) return(NULL)
  mcp_probesets <- load_data("mcp_probesets")
  if (is.null(mcp_probesets)) return(NULL)

  res <- MCPcounter.estimate(
    eset,
    featuresType = "HUGO_symbols",
    probesets = mcp_probesets,
    genes = mcp_genes
  )
  res <- as.data.frame(t(res))

  .format_deconv_result(res, project, "MCPcounter")
}

#' Deconvolve Immune Microenvironment Using EPIC
#'
#' @description
#' Estimates immune cell fractions using EPIC algorithm.
#'
#' @param eset Gene expression matrix with genes as row names.
#' @param project Optional project name. Default is `NULL`.
#' @param tumor Logical indicating tumor (`TRUE`) or normal (`FALSE`) samples.
#'
#' @return Data frame with EPIC cell fraction estimates. Columns suffixed with
#'   `_EPIC`.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' TRef <- load_data("TRef")
#' if (!is.null(TRef)) {
#'   set.seed(123)
#'   sim_eset <- matrix(rnorm(nrow(TRef$refProfiles) * 3), nrow(TRef$refProfiles), 3)
#'   rownames(sim_eset) <- rownames(TRef$refProfiles)
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'   
#'   # Run deconvolution
#'   result <- deconvo_epic(eset = sim_eset, project = "Example", tumor = TRUE)
#'   if (!is.null(result)) head(result)
#' }
deconvo_epic <- function(eset, project = NULL, tumor = TRUE) {
  if (is.null(eset)) return(NULL)
  cli::cli_alert_info("Running EPIC deconvolution")

  ref <- if (tumor) "TRef" else "BRef"

  out <- EPIC(
    bulk = eset,
    reference = ref,
    mRNA_cell = NULL,
    scaleExprs = TRUE
  )
  if (is.null(out)) return(NULL)
  res <- as.data.frame(out$cellFractions)

  .format_deconv_result(res, project, "EPIC")
}

#' Deconvolve Using CIBERSORT
#'
#' @description
#' CIBERSORT is freely available to academic users. License and binary can be
#' obtained from https://cibersortx.stanford.edu.
#'
#' @param eset Expression matrix with gene symbols as row names.
#' @param project Optional project name. Default is `NULL`.
#' @param arrays Logical: optimized for microarray data. Default is `FALSE`.
#' @param perm Permutations for statistical analysis. Default is 1000.
#' @param absolute Logical: run in absolute mode. Default is `FALSE`.
#' @param abs_method Method for absolute mode: `"sig.score"` or `"no.sumto1"`.
#'   Default is `"sig.score"`.
#' @param parallel Enable parallel execution. Default is `FALSE`.
#' @param num_cores Number of cores for parallel mode. Default is 2.
#' @param seed Random seed for reproducibility. Default is `NULL`.
#'
#' @return Data frame with CIBERSORT cell fractions. Columns suffixed with `_CIBERSORT`.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' lm22 <- load_data("lm22")
#' if (!is.null(lm22)) {
#'   set.seed(123)
#'   sim_eset <- matrix(rnorm(nrow(lm22) * 3), nrow(lm22), 3)
#'   rownames(sim_eset) <- rownames(lm22)
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'
#'   # Run deconvolution
#'   result <- deconvo_cibersort(
#'     eset = sim_eset,
#'     project = "TCGA-STAD",
#'     perm = 10
#'   )
#'   if (!is.null(result)) head(result)
#' }
deconvo_cibersort <- function(eset,
                              project = NULL,
                              arrays = FALSE,
                              perm = 1000,
                              absolute = FALSE,
                              abs_method = "sig.score",
                              parallel = FALSE,
                              num_cores = 2,
                              seed = NULL) {
  if (is.null(eset)) return(NULL)
  mode_label <- if (absolute) "CIBERSORT (absolute mode)" else "CIBERSORT"
  cli::cli_alert_info("Running {mode_label}")

  eset <- as.data.frame(eset)
  quantile_norm <- arrays # Disable for RNA-seq

  sig_matrix <- load_data("lm22")
  if (is.null(sig_matrix)) return(NULL)

  res <- CIBERSORT(
    sig_matrix = sig_matrix,
    mixture_file = eset,
    perm = perm,
    QN = quantile_norm,
    absolute = absolute,
    abs_method = abs_method,
    parallel = parallel,
    num_cores = num_cores,
    seed = seed
  )
  if (is.null(res)) return(NULL)

  .format_deconv_result(as.data.frame(res), project, "CIBERSORT")
}

#' Calculate Immunophenoscore (IPS)
#'
#' @description
#' Calculates immune phenotype scores from gene expression data.
#'
#' @param eset Gene expression matrix.
#' @param project Optional project name. Default is `NULL`.
#' @param plot Logical: generate visualization. Default is `FALSE`.
#'
#' @return Data frame with IPS scores. Columns suffixed with `_IPS`.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' ips_genes <- load_data("ips_gene_set")
#' if (!is.null(ips_genes)) {
#'   set.seed(123)
#'   sim_eset <- matrix(rnorm(nrow(ips_genes) * 3), nrow(ips_genes), 3)
#'   rownames(sim_eset) <- ips_genes$GENE
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'   
#'   # Run calculation
#'   result <- deconvo_ips(eset = sim_eset, project = "Example")
#'   if (!is.null(result)) head(result)
#' }
deconvo_ips <- function(eset, project = NULL, plot = FALSE) {
  if (is.null(eset)) return(NULL)
  cli::cli_alert_info("Running IPS calculation")

  res <- IPS_calculation(project = project, eset = eset, plot = plot)
  if (is.null(res)) return(NULL)

  .format_deconv_result(res, project, "IPS")
}

#' Calculate ESTIMATE Scores
#'
#' @description
#' Calculates stromal, immune, and ESTIMATE scores from gene expression.
#'
#' @param eset Gene expression matrix with gene symbols.
#' @param project Optional project name. Default is `NULL`.
#' @param platform Platform type: `"affymetrix"` or `"illumina"`.
#'   Default is `"affymetrix"`.
#'
#' @return Data frame with ESTIMATE scores. Columns suffixed with `_estimate`.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' common_genes <- load_data("common_genes")
#' if (!is.null(common_genes)) {
#'   set.seed(123)
#'   sim_eset <- matrix(rnorm(nrow(common_genes) * 3), nrow(common_genes), 3)
#'   rownames(sim_eset) <- common_genes$GeneSymbol
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'   
#'   # Run calculation
#'   result <- deconvo_estimate(sim_eset, project = "TCGA-STAD")
#'   if (!is.null(result)) head(result)
#' }
deconvo_estimate <- function(eset, project = NULL, platform = "affymetrix") {
  if (is.null(eset)) return(NULL)
  cli::cli_alert_info("Running ESTIMATE")

  # Save original sample IDs for later restoration
  original_sample_ids <- colnames(eset)

  eset <- as.data.frame(eset)
  eset <- tibble::rownames_to_column(eset, var = "symbol")

  sampleData <- tempfile(pattern = "estimate_", fileext = ".txt")
  on.exit(unlink(sampleData), add = TRUE)

  utils::write.table(eset, sampleData,
    sep = "\t",
    row.names = FALSE, quote = FALSE
  )

  gct_file <- tempfile(pattern = "estimate_", fileext = ".gct")
  score_file <- tempfile(pattern = "estimate_", fileext = ".gct")
  on.exit(unlink(c(gct_file, score_file)), add = TRUE)

  filter_res <- filterCommonGenes(
    input.f = sampleData,
    output.f = gct_file,
    id = "GeneSymbol"
  )
  if (isFALSE(filter_res)) return(NULL)

  score_res <- estimateScore(
    input.ds = gct_file,
    output.ds = score_file,
    platform = platform
  )
  if (isFALSE(score_res)) return(NULL)

  if (!file.exists(score_file)) return(NULL)
  scores <- tryCatch(utils::read.table(score_file,
    skip = 2, header = TRUE,
    stringsAsFactors = FALSE
  ), error = function(e) NULL)
  if (is.null(scores)) return(NULL)
  
  rownames(scores) <- scores[, 1]
  scores <- t(scores[, 3:ncol(scores), drop = FALSE])

  # Restore original sample IDs (R's read.table converts - to . in column names)
  rownames(scores) <- original_sample_ids

  .format_deconv_result(as.data.frame(scores), project, "estimate")
}

#' Deconvolve Using Custom Reference
#'
#' @description
#' Cell fraction estimation using SVR or lsei methods with custom reference.
#'
#' @param eset Gene expression matrix.
#' @param project Optional project name. Default is `NULL`.
#' @param arrays Logical: use quantile normalization. Default is `TRUE`.
#' @param method Method: `"svr"` or `"lsei"`. Default is `"svr"`.
#' @param perm Permutations for SVR. Default is 100.
#' @param reference Custom reference matrix (e.g., lm22, lm6).
#' @param scale_reference Logical: scale reference. Default is `TRUE`.
#' @param absolute.mode Logical: absolute mode for SVR. Default is `FALSE`.
#' @param abs.method Method for absolute mode. Default is `"sig.score"`.
#'
#' @return Data frame with cell fractions. Columns suffixed with `_CIBERSORT`.
#'
#' @author Dongqiang Zeng, Rongfang Shen
#' @export
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' sim_ref <- matrix(rnorm(100 * 5), 100, 5)
#' rownames(sim_ref) <- paste0("Gene", 1:100)
#' colnames(sim_ref) <- paste0("CellType", 1:5)
#'
#' sim_eset <- matrix(rnorm(100 * 3), 100, 3)
#' rownames(sim_eset) <- paste0("Gene", 1:100)
#' colnames(sim_eset) <- paste0("Sample", 1:3)
#'
#' # Run deconvolution
#' result <- deconvo_ref(eset = sim_eset, reference = sim_ref, method = "lsei")
#' if (!is.null(result)) head(result)
deconvo_ref <- function(eset,
                        project = NULL,
                        arrays = TRUE,
                        method = c("svr", "lsei"),
                        perm = 100,
                        reference,
                        scale_reference = TRUE,
                        absolute.mode = FALSE,
                        abs.method = "sig.score") {
  if (is.null(eset)) return(NULL)
  method <- rlang::arg_match(method)
  rlang::check_installed("limSolve")

  # Check gene overlap
  common_genes <- intersect(rownames(eset), rownames(reference))
  if (length(common_genes) == 0) {
    cli::cli_abort(c(
      "No matching genes between eset and reference.",
      "i" = "Check gene identifier formats match."
    ))
  }
  cli::cli_alert_info("Found {length(common_genes)} common genes")

  if (method == "svr") {
    cli::cli_alert_info("Running SVR deconvolution")

    res <- CIBERSORT(
      sig_matrix = reference,
      mixture_file = as.data.frame(eset),
      perm = perm,
      QN = arrays,
      absolute = absolute.mode,
      abs_method = abs.method
    )
    res <- as.data.frame(res)
  } else if (method == "lsei") {
    cli::cli_alert_info("Running lsei deconvolution")

    eset <- as.matrix(eset)
    reference <- as.matrix(reference)

    # Quantile normalization
    if (arrays) {
      rlang::check_installed("preprocessCore")
      eset <- preprocessCore::normalize.quantiles(eset)
    }

    # Scale reference (global scaling for backward compatibility)
    if (scale_reference) {
      reference <- (reference - mean(reference)) / sd(as.vector(reference))
    }

    # Find common genes
    common <- intersect(rownames(eset), rownames(reference))
    eset <- eset[match(common, rownames(eset)), , drop = FALSE]
    reference <- reference[match(common, rownames(reference)), , drop = FALSE]

    # Setup lsei constraints
    Numofx <- ncol(reference)
    AA <- reference
    EE <- rep(1, Numofx)
    FF <- 1
    GG <- diag(nrow = Numofx)
    HH <- rep(0, Numofx)

    # Run deconvolution
    out.all <- vapply(seq_len(ncol(eset)), function(i) {
      BB <- eset[, i]
      BB <- (BB - mean(BB)) / stats::sd(BB)
      out <- limSolve::lsei(AA, BB, EE, FF, GG, HH)
      out$X
    }, numeric(Numofx))

    res <- as.data.frame(t(out.all))
    colnames(res) <- colnames(reference)
    rownames(res) <- colnames(eset)
  }

  .format_deconv_result(res, project, "CIBERSORT")
}

#' Deconvolve Using TIMER
#'
#' @description
#' TIMER deconvolution for cancer-specific immune estimation.
#'
#' @param eset Gene expression matrix.
#' @param project Optional project name. Default is `NULL`.
#' @param indications Cancer type for each sample (e.g., `"brca"`, `"stad"`).
#'   Must match number of columns in `eset`.
#'
#' @return Data frame with TIMER cell fractions. Columns suffixed with `_TIMER`.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' immune <- load_data("immuneCuratedData")
#' cancer_genes <- load_data("cancer_type_genes")
#' if (!is.null(immune) && !is.null(cancer_genes)) {
#'   set.seed(123)
#'   # Extract a subset of genes to avoid "NA" gene names breaking read.table
#'   genes <- unique(c(head(rownames(immune$genes), 100), cancer_genes[["stad"]]))
#'   sim_eset <- matrix(rnorm(length(genes) * 3), length(genes), 3)
#'   rownames(sim_eset) <- genes
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'   
#'   # Run deconvolution
#'   result <- deconvo_timer(
#'     eset = sim_eset, project = "TCGA-STAD",
#'     indications = rep("stad", 3)
#'   )
#'   if (!is.null(result)) head(result)
#' }
deconvo_timer <- function(eset, project = NULL, indications = NULL) {
  if (is.null(eset)) return(NULL)
  cli::cli_alert_info("Running TIMER deconvolution")

  if (!is.null(indications)) {
    indications <- tolower(indications)

    if (length(indications) != ncol(eset)) {
      cli::cli_abort(c(
        "Length of 'indications' must match number of samples.",
        "i" = "Got {length(indications)} indications for {ncol(eset)} samples"
      ))
    }
  }

  # Prepare temp files
  args <- new.env()
  args$outdir <- tempdir()
  args$batch <- tempfile()

  # Write data for each cancer type
  unique_inds <- unique(indications)
  for (ind in unique_inds) {
    tmp_file <- tempfile()
    tmp_mat <- eset[, indications == ind, drop = FALSE]
    tmp_mat <- tibble::as_tibble(tmp_mat, rownames = "gene_symbol")
    utils::write.table(tmp_mat, tmp_file,
      sep = "\t",
      quote = FALSE, row.names = FALSE
    )
    cat(paste0(tmp_file, ",", ind, "\n"), file = args$batch, append = TRUE)
  }

  # Run TIMER
  results <- deconvolute_timer.default(args)
  if (is.null(results)) return(NULL)
  results <- results[, colnames(eset), drop = FALSE]
  colnames(results) <- colnames(eset)
  results <- as.data.frame(t(results))

  .format_deconv_result(results, project, "TIMER")
}

#' Deconvolve Using quanTIseq
#'
#' @description
#' quanTIseq deconvolution for RNA-seq immune cell fractions.
#'
#' @param eset Gene expression matrix.
#' @param project Optional project name. Default is `NULL`.
#' @param tumor Logical: tumor samples. Must be specified.
#' @param arrays Logical: microarray data. Must be specified.
#' @param scale_mrna Logical: correct for mRNA content. Must be specified.
#'
#' @return Data frame with quanTIseq cell fractions. Columns suffixed with
#'   `_quantiseq`.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' quantiseq_data <- load_data("quantiseq_data")
#' if (!is.null(quantiseq_data)) {
#'   set.seed(123)
#'   n_sig <- nrow(quantiseq_data$TIL10_signature)
#'   sim_eset <- matrix(rnorm(n_sig * 3), n_sig, 3)
#'   rownames(sim_eset) <- rownames(quantiseq_data$TIL10_signature)
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'   
#'   # Run deconvolution
#'   result <- deconvo_quantiseq(eset = sim_eset, project = "Example", tumor = TRUE,
#'                               arrays = FALSE, scale_mrna = FALSE)
#'   if (!is.null(result)) head(result)
#' }
deconvo_quantiseq <- function(eset, project = NULL, tumor, arrays, scale_mrna) {
  if (is.null(eset)) return(NULL)
  cli::cli_alert_info("Running quanTIseq deconvolution")

  res <- deconvolute_quantiseq.default(
    mix.mat = eset,
    tumor = tumor,
    arrays = arrays,
    mRNAscale = scale_mrna
  )
  if (is.null(res)) return(NULL)

  res <- tibble::column_to_rownames(tibble::as_tibble(res), var = "Sample")

  .format_deconv_result(res, project, "quantiseq")
}

#' Main TME Deconvolution Function
#'
#' @description
#' Unified interface for multiple TME deconvolution methods.
#'
#' @param eset Gene expression matrix with HGNC symbols as row names.
#' @param project Optional project name. Default is `NULL`.
#' @param method Deconvolution method. See [tme_deconvolution_methods].
#' @param arrays Logical: microarray-optimized mode. Default is `FALSE`.
#' @param tumor Logical: tumor-optimized mode (EPIC). Default is `TRUE`.
#' @param perm Permutations (CIBERSORT/SVR). Default is 1000.
#' @param reference Custom reference matrix (SVR/lsei).
#' @param scale_reference Logical: scale reference (SVR/lsei).
#' @param plot Logical: generate plots (IPS). Default is `FALSE`.
#' @param scale_mrna Logical: mRNA correction (quanTIseq/EPIC).
#' @param group_list Cancer types for TIMER (vector).
#' @param platform Platform for ESTIMATE. Default is `"affymetrix"`.
#' @param absolute.mode Logical: absolute mode (CIBERSORT/SVR).
#'   Default is `FALSE`.
#' @param abs.method Absolute mode method. Default is `"sig.score"`.
#' @param ... Additional arguments passed to method.
#'
#' @return Tibble with cell fractions and `ID` column.
#'
#' @author Dongqiang Zeng, Rongfang Shen
#' @export
#'
#' @references
#' \enumerate{
#'   \item Newman et al. (2015). Robust enumeration of cell subsets from tissue
#'     expression profiles. Nature Methods.
#'   \item Vegesna et al. (2013). Inferring tumour purity and stromal/immune
#'     cell admixture. Nature Communications.
#'   \item Finotello et al. (2019). Molecular and pharmacological modulators of
#'     the tumor immune contexture. Genome Medicine.
#'   \item Li et al. (2016). Comprehensive analyses of tumor immunity.
#'     Genome Biology.
#'   \item Charoentong et al. (2017). Pan-cancer Immunogenomic Analyses.
#'     Cell Reports.
#'   \item Becht et al. (2016). Estimating population abundance of tissue-infiltrating
#'     immune cells. Genome Biology.
#'   \item Aran et al. (2017). xCell: digitally portraying tissue cellular
#'     heterogeneity. Genome Biology.
#'   \item Racle et al. (2017). Simultaneous enumeration of cancer and immune
#'     cell types. ELife.
#' }
#'
#' @examples
#' mcp_genes <- load_data("mcp_genes")
#' if (!is.null(mcp_genes)) {
#'   set.seed(123)
#'   sim_eset <- matrix(rnorm(nrow(mcp_genes) * 3), nrow(mcp_genes), 3)
#'   rownames(sim_eset) <- mcp_genes$`HUGO symbols`
#'   colnames(sim_eset) <- paste0("Sample", 1:3)
#'   
#'   # Run deconvolution
#'   result <- deconvo_tme(eset = sim_eset, method = "mcpcounter")
#'   if (!is.null(result)) head(result)
#' }
deconvo_tme <- function(eset,
                        project = NULL,
                        method = tme_deconvolution_methods,
                        arrays = FALSE,
                        tumor = TRUE,
                        perm = 1000,
                        reference,
                        scale_reference = TRUE,
                        plot = FALSE,
                        scale_mrna = TRUE,
                        group_list = NULL,
                        platform = "affymetrix",
                        absolute.mode = FALSE,
                        abs.method = "sig.score",
                        ...) {
  if (is.null(eset)) return(NULL)
  method <- rlang::arg_match(method, tme_deconvolution_methods)

  # Validate input
  if (any(grepl("^ENSG000", rownames(eset)))) {
    cli::cli_abort(c(
      "Ensembl IDs detected in row names.",
      "i" = "Most deconvolution methods require HGNC gene symbols.",
      "*" = "Use anno_eset() to convert Ensembl IDs to gene symbols"
    ))
  }

  if (max(eset, na.rm = TRUE) < 50) {
    cli::cli_warn(c(
      "Data values appear small (< 50).",
      "i" = "Input should be in TPM/FPKM scale, not log-transformed"
    ))
  }

  # Dispatch to method
  res <- switch(method,
    xcell = deconvo_xcell(eset, project, arrays = arrays, ...),
    mcpcounter = deconvo_mcpcounter(eset, project, ...),
    epic = deconvo_epic(eset, project, tumor = tumor, ...),
    cibersort = deconvo_cibersort(eset, project,
      absolute = absolute.mode,
      arrays = arrays, perm = perm, ...
    ),
    cibersort_abs = deconvo_cibersort(eset, project,
      absolute = TRUE,
      abs_method = abs.method,
      arrays = arrays, perm = perm, ...
    ),
    ips = deconvo_ips(eset, project, plot = plot, ...),
    quantiseq = deconvo_quantiseq(eset, project,
      tumor = tumor,
      arrays = arrays, scale_mrna = scale_mrna, ...
    ),
    estimate = deconvo_estimate(eset, project, platform, ...),
    timer = deconvo_timer(eset, project, indications = group_list, ...),
    svr = deconvo_ref(eset, project,
      reference = reference,
      arrays = arrays, method = "svr",
      absolute.mode = absolute.mode, ...
    ),
    lsei = deconvo_ref(eset, project,
      reference = reference,
      arrays = arrays, method = "lsei",
      scale_reference = scale_reference, ...
    )
  )

  if (is.null(res)) return(NULL)
  tibble::as_tibble(res)
}
