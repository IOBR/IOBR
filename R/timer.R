#' Source code for the TIMER deconvolution method.
#'
#' This code is adapted from https://github.com/hanfeisun/TIMER, which
#' again is an adapted version of the original TIMER source code
#' from http://cistrome.org/TIMER/download.html.
#'
#' The method is described in Li et al. Genome Biology 2016;17(1):174., PMID
#' 27549193.


#' Display Timer Information Messages
#'
#' @description
#' Formats and displays informational messages for timing or logging purposes.
#' Useful for tracking progress or stages of execution within scripts.
#'
#' @param string Character. Message to be displayed.
#'
#' @return None; used for its side effect of printing a message.
#'
#' @export
#' @author Bo Li
#'
#' @examples
#' timer_info("Data processing started.")
timer_info <- function(string) {
  cli::cli_alert_info(string)
}


#' TIMER Available Cancer Types
#'
#' Character vector of cancer types supported by TIMER deconvolution.
#' TIMER signatures are cancer-specific.
#'
#' @return Character vector of available cancer type abbreviations.
#'
#' @export
#'
#' @examples
#' # List all available cancer types for TIMER
#' timer_available_cancers
#'
#' # Check if a cancer type is supported
#' "brca" %in% timer_available_cancers
timer_available_cancers <- c(
  "kich", "blca", "brca", "cesc", "gbm", "hnsc", "kirp", "lgg",
  "lihc", "luad", "lusc", "prad", "sarc", "pcpg", "paad", "tgct",
  "ucec", "ov", "skcm", "dlbc", "kirc", "acc", "meso", "thca",
  "uvm", "ucs", "thym", "esca", "stad", "read", "coad", "chol"
)


#' Remove Batch Effect of Expression Set
#'
#' @description
#' Removes batch effects between two gene expression datasets, typically
#' representing different sample types such as cancer cells and immune cells.
#' Uses ComBat from the sva package for batch correction.
#'
#' @param cancer.exp Matrix or data frame. Cancer cell expression data with
#'   genes as rows and samples as columns.
#' @param immune.exp Matrix or data frame. Immune cell expression data with
#'   genes as rows and samples as columns.
#' @param immune.cellType Vector. Cell type for each column in `immune.exp`.
#'
#' @return A list containing:
#' \describe{
#'   \item{1}{Batch effect corrected cancer expression data}
#'   \item{2}{Batch effect corrected immune expression data}
#'   \item{3}{Aggregated immune expression data (median per cell type)}
#' }
#'
#' @export
#' @author Bo Li
#'
#' @examples
#' set.seed(123)
#' gene_names <- paste0("Gene", 1:100)
#' sample_names_cancer <- paste0("CancerSample", 1:10)
#' cancer.exp <- matrix(runif(1000, 1, 1000),
#'   nrow = 100, ncol = 10,
#'   dimnames = list(gene_names, sample_names_cancer)
#' )
#'
#' sample_names_immune <- paste0("ImmuneSample", 1:5)
#' immune.exp <- matrix(runif(500, 1, 1000),
#'   nrow = 100, ncol = 5,
#'   dimnames = list(gene_names, sample_names_immune)
#' )
#'
#' immune.cellType <- c("T-cell", "B-cell", "T-cell", "NK-cell", "B-cell")
#' names(immune.cellType) <- sample_names_immune
#' \donttest{
#' result <- RemoveBatchEffect(cancer.exp, immune.exp, immune.cellType)
#' }
RemoveBatchEffect <- function(cancer.exp, immune.exp, immune.cellType) {
  rlang::check_installed("sva")

  tmp.dd <- as.matrix(cancer.exp)
  tmp.ss <- intersect(rownames(tmp.dd), rownames(immune.exp))

  N1 <- ncol(tmp.dd)
  tmp.dd <- cbind(tmp.dd[tmp.ss, ], immune.exp[tmp.ss, ])
  tmp.dd <- as.matrix(tmp.dd)
  mode(tmp.dd) <- "numeric"

  N2 <- ncol(immune.exp)
  tmp.batch <- c(rep(1, N1), rep(2, N2))

  rlang::check_installed(c("sva", "BiocParallel"))
  suppressMessages(
    tmp.dd0 <- sva::ComBat(tmp.dd, tmp.batch, c(),
      BPPARAM = BiocParallel::bpparam("SerialParam")
    )
  )

  dd.br <- tmp.dd0[, 1:N1]
  immune.exp.br <- tmp.dd0[, (N1 + 1):(N1 + N2)]

  tmp0 <- NULL
  for (kk in unique(names(immune.cellType))) {
    tmp.vv <- which(names(immune.cellType) == kk)

    if (length(tmp.vv) == 1) {
      median_expression <- apply(immune.exp.br[, tmp.vv, drop = FALSE], 1, median, na.rm = TRUE)
    } else {
      median_expression <- apply(immune.exp.br[, tmp.vv], 1, median, na.rm = TRUE)
    }

    if (is.null(tmp0)) {
      tmp0 <- median_expression
    } else {
      tmp0 <- cbind(tmp0, median_expression)
    }
  }

  immune.exp.agg.br <- tmp0
  colnames(immune.exp.agg.br) <- unique(names(immune.cellType))

  list(as.matrix(dd.br), immune.exp.br, immune.exp.agg.br)
}


#' Process Batch Table and Validate Cancer Types
#'
#' @description
#' Processes input data containing cancer types and validates each category
#' against a predefined list of supported cancer types (`timer_available_cancers`).
#'
#' @param args A list containing input parameters:
#' \describe{
#'   \item{batch}{Character. Path to a CSV file (optional).}
#'   \item{expression}{Character vector of expression identifiers (used if
#'     `batch` is NULL).}
#'   \item{category}{Character vector of cancer types (used if `batch` is NULL).}
#' }
#'
#' @return A character matrix with two columns:
#' \describe{
#'   \item{Column 1}{Expression identifiers}
#'   \item{Column 2}{Cancer categories}
#' }
#'
#' @export
#'
#' @examples
#' args <- list(
#'   expression = c("exp1", "exp2"),
#'   category = c("luad", "brca"),
#'   batch = NULL
#' )
#' result <- check_cancer_types(args)
check_cancer_types <- function(args) {
  if (!is.null(args$batch)) {
    timer_info("Enter batch mode")
    cancers <- as.matrix(read.table(args$batch, sep = ",", stringsAsFactors = FALSE))
  } else {
    if (length(args$expression) != length(args$category)) {
      cli::cli_abort("expression and category must have the same length")
    }
    cancers <- cbind(args$expression, args$category)
  }

  invalid_cancers <- setdiff(cancers[, 2], timer_available_cancers)
  if (length(invalid_cancers) > 0) {
    cli::cli_abort("Unknown cancer type(s): {.val {invalid_cancers}}")
  }

  cancers
}


#' Constrained Regression Method (Abbas et al., 2009)
#'
#' @description
#' Implements a constrained regression approach described by Abbas et al. (2009).
#' Estimates proportions of immune cell types within mixed cancer tissue samples
#' based on gene expression data. Iteratively adjusts regression coefficients
#' to ensure non-negative values.
#'
#' @param XX Matrix. Immune expression data with genes as rows and cell types
#'   as columns.
#' @param YY Vector. Cancer expression data with gene expression levels.
#' @param w Vector or NA. Weights for regression. Default is NA (no weights).
#'
#' @return Vector with non-negative coefficients representing proportions of
#'   each cell type.
#'
#' @export
#'
#' @examples
#' XX <- matrix(runif(100), nrow = 10, ncol = 10)
#' colnames(XX) <- paste("CellType", 1:10, sep = "")
#' YY <- runif(10)
#' results <- GetFractions.Abbas(XX, YY)
#' print(results)
GetFractions.Abbas <- function(XX, YY, w = NA) {
  ss.remove <- c()
  ss.names <- colnames(XX)

  while (TRUE) {
    if (length(ss.remove) == 0) {
      tmp.XX <- XX
    } else {
      if (is.null(ncol(tmp.XX))) {
        return(rep(0, ncol(XX)))
      }
      tmp.XX <- tmp.XX[, -ss.remove]
    }

    if (length(ss.remove) > 0) {
      ss.names <- ss.names[-ss.remove]
      if (length(ss.names) == 0) {
        return(rep(0, ncol(XX)))
      }
    }

    if (is.na(w[1])) {
      tmp <- lsfit(tmp.XX, YY, intercept = FALSE)
    } else {
      tmp <- lsfit(tmp.XX, YY, w, intercept = FALSE)
    }

    if (is.null(ncol(tmp.XX))) {
      tmp.beta <- tmp$coefficients[1]
    } else {
      tmp.beta <- tmp$coefficients[seq_len(ncol(tmp.XX))]
    }

    if (min(tmp.beta) > 0) break
    ss.remove <- which.min(tmp.beta)
  }

  tmp.F <- rep(0, ncol(XX))
  names(tmp.F) <- colnames(XX)
  tmp.F[ss.names] <- tmp.beta
  tmp.F
}

#' Convert Rowname To Loci
#'
#' @description
#' Processes a gene expression data matrix by modifying its row names.
#' Extracts the gene identifier from row names formatted as 'GENE|ID',
#' simplifying them to 'GENE'.
#'
#' @param cancerGeneExpression Matrix or data frame. Gene expression data with
#'   row names in the format 'GENE|ID'.
#'
#' @return Matrix with modified gene expression data with updated row names.
#'   Rows without a valid identifier are removed.
#'
#' @export
#'
#' @examples
#' example_data <- matrix(runif(20), ncol = 5)
#' rownames(example_data) <- c("LOC101", "LOC102", "LOC103", "LOC104")
#' processed_data <- ConvertRownameToLoci(example_data)
#' print(processed_data)
ConvertRownameToLoci <- function(cancerGeneExpression) {
  tmp <- strsplit(rownames(cancerGeneExpression), "\\|")
  tmp <- sapply(tmp, function(x) x[[1]])
  tmp.vv <- which(nchar(tmp) > 1)
  rownames(cancerGeneExpression) <- tmp
  extracted <- as.matrix(cancerGeneExpression[tmp.vv, ])
  colnames(extracted) <- colnames(cancerGeneExpression)
  extracted
}


#' Parse Input Gene Expression Data
#'
#' @description
#' Reads gene expression data from a tab-delimited text file, using the first
#' column as row names. Converts data into a numeric matrix for analysis.
#'
#' @param path Character. Path to a tab-delimited gene expression file.
#'   First column should contain gene identifiers.
#'
#' @return Numeric matrix of gene expression values with genes as rows and
#'   samples as columns.
#'
#' @export
#'
#' @examples
#' tf <- tempfile(fileext = ".tsv")
#' expr <- data.frame(
#'   gene = c("GeneA", "GeneB", "GeneC"),
#'   Sample1 = c(10, 20, 30),
#'   Sample2 = c(15, 25, 35)
#' )
#' write.table(expr, tf, sep = "\t", row.names = FALSE, quote = FALSE)
#' gene_expression_data <- ParseInputExpression(tf)
#' print(gene_expression_data)
ParseInputExpression <- function(path) {
  ret <- read.delim(path, row.names = 1, check.names = FALSE)
  ret <- as.matrix(ret)
  mode(ret) <- "numeric"
  ret
}


#' Draw QQ Plot Comparing Cancer and Immune Expression
#'
#' @description
#' Creates a quantile-quantile (QQ) plot to compare gene expression
#' distributions between cancer and immune samples. Points along the diagonal
#' indicate similar distributions.
#'
#' @param cancer.exp Vector. Gene expression data for cancer samples.
#' @param immune.exp Vector. Gene expression data for immune samples.
#' @param name Character. Optional subtitle with additional information.
#'
#' @return Generates a QQ plot.
#'
#' @export
#'
#' @examples
#' cancer_exp <- rnorm(100, mean = 5, sd = 1.5)
#' immune_exp <- rnorm(100, mean = 5, sd = 1.5)
#' DrawQQPlot(
#'   cancer.exp = cancer_exp,
#'   immune.exp = immune_exp,
#'   name = "Comparison of Gene Expression"
#' )
DrawQQPlot <- function(cancer.exp, immune.exp, name = "") {
  qq <- qqplot(cancer.exp, immune.exp,
    xlab = "Tumor Expression", ylab = "Ref Expression",
    main = "Sample-Sample Q-Q plot"
  )
  mtext(name, col = "gray11")

  start <- ceiling(0.4 * length(qq$x))
  end <- floor(0.9 * length(qq$x))
  qq.sub <- list(x = qq$x[start:end], y = qq$y[start:end])
  fit <- lm(y ~ x, data = qq.sub)
  abline(fit, col = "blue")
  invisible(qq)
}

#' Get Outlier Genes
#'
#' @description
#' Identifies outlier genes from multiple cancer datasets. Treats the top 5
#' expressed genes in each sample as outliers and returns unique outlier genes.
#'
#' @param cancers Data frame. One column containing paths to gene expression files.
#'
#' @return Vector of unique gene names identified as outliers.
#'
#' @export
#'
#' @examples
#' tf1 <- tempfile(fileext = ".tsv")
#' tf2 <- tempfile(fileext = ".tsv")
#'
#' expr1 <- data.frame(
#'   gene = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"),
#'   Sample1 = c(10, 50, 30, 80, 60, 20),
#'   Sample2 = c(15, 40, 25, 90, 55, 10)
#' )
#'
#' expr2 <- data.frame(
#'   gene = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"),
#'   Sample3 = c(100, 20, 10, 60, 30, 80),
#'   Sample4 = c(95, 25, 15, 70, 35, 85)
#' )
#'
#' write.table(expr1, tf1, sep = "\t", row.names = FALSE, quote = FALSE)
#' write.table(expr2, tf2, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' cancers <- data.frame(ExpressionFiles = c(tf1, tf2))
#' outlier_genes <- GetOutlierGenes(cancers)
GetOutlierGenes <- function(cancers) {
  outlier.total <- c()

  for (i in seq_len(nrow(cancers))) {
    cancer.expFile <- as.character(cancers[i, 1])
    cancer.expression <- ParseInputExpression(cancer.expFile)

    for (j in seq_len(ncol(cancer.expression))) {
      outlier <- rownames(cancer.expression)[tail(order(cancer.expression[, j]), 5)]
      outlier.total <- c(outlier.total, outlier)
    }
  }

  unique(outlier.total)
}


#' Deconvolute Tumor Microenvironment Using TIMER
#'
#' @description
#' Performs deconvolution of the tumor microenvironment using the TIMER
#' algorithm. Processes multiple cancer datasets, removes batch effects,
#' and estimates immune cell type abundances.
#'
#' @param args List or environment containing parameters:
#' \describe{
#'   \item{outdir}{Character. Output directory path.}
#'   \item{batch}{Character. File containing paths to expression data and
#'     cancer types.}
#' }
#'
#' @return Matrix of abundance scores for different immune cell types across
#'   multiple cancer samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # file
#' tf <- tempfile(fileext = ".csv")
#' write.table(data.frame("exp1", "luad", "exp2", "brca"),
#'   file = tf, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
#' )
#' outdir <- tempdir()
#' args <- list(outdir = outdir, batch = tf)
#' results <- deconvolute_timer.default(args)
#' }
deconvolute_timer.default <- function(args) {
  cancers <- check_cancer_types(args)

  timer_info("Loading immune gene expression")

  immune <- load_data("immuneCuratedData")
  immune.geneExpression <- immune$genes
  immune.cellTypes <- immune$celltypes

  outlier.genes <- sort(GetOutlierGenes(cancers))
  cli::cli_alert_info("Outlier genes: {paste(outlier.genes, collapse = ' ')}")

  dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
  results_dir <- file.path(args$outdir, "results")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }

  abundance.score.matrix <- c()
  pdf(file.path(results_dir, "output.pdf"))

  for (i in seq_len(nrow(cancers))) {
    cancer.expFile <- cancers[i, 1]
    cancer.category <- cancers[i, 2]

    cancer.expression <- ParseInputExpression(cancer.expFile)
    index <- !(row.names(cancer.expression) %in% outlier.genes)
    cancer.expression <- cancer.expression[index, , drop = FALSE]
    cancer.colnames <- colnames(cancer.expression)

    timer_info(paste("Removing batch effects for", cancer.category))

    for (j in seq_along(cancer.colnames)) {
      DrawQQPlot(cancer.expression[, j], immune.geneExpression[, 1], name = cancer.colnames[j])
    }

    tmp <- RemoveBatchEffect(cancer.expression, immune.geneExpression, immune.cellTypes)
    cancer.expNorm <- tmp[[1]]
    immune.expNormMedian <- tmp[[3]]

    for (j in seq_along(cancer.colnames)) {
      DrawQQPlot(cancer.expNorm[, j], immune.expNormMedian[, 1],
        name = paste("After batch removing and aggregating for", cancer.colnames[j])
      )
    }

    cancer_type_genes_data <- load_data("cancer_type_genes")
    gene.selected.marker <- cancer_type_genes_data[[which(names(cancer_type_genes_data) == cancer.category)]]
    gene.selected.marker <- intersect(gene.selected.marker, row.names(cancer.expNorm))

    if (length(gene.selected.marker) < 6) {
      cli::cli_abort(c(
        "Insufficient marker genes for TIMER deconvolution.",
        "i" = "Cancer type '{cancer.category}' requires at least 6 marker genes.",
        "i" = "Found only {length(gene.selected.marker)} marker genes in the input data.",
        "*" = "Use the full expression matrix instead of a subset (e.g., eset[1:500, ]).",
        "*" = "TIMER requires cancer-specific gene markers that may not be in top-expressed genes."
      ))
    }

    XX <- immune.expNormMedian[gene.selected.marker, -4]
    YY <- cancer.expNorm[gene.selected.marker, , drop = FALSE]

    for (j in seq_along(cancer.colnames)) {
      fractions <- GetFractions.Abbas(XX, YY[, j])
      barplot(fractions,
        cex.names = 0.8, names.arg = names(fractions),
        xlab = "cell type", ylab = "abundance",
        main = paste("Abundance estimation for", cancer.colnames[j])
      )
      box()

      abundance.score.matrix <- cbind(abundance.score.matrix, fractions)
      colnames(abundance.score.matrix)[ncol(abundance.score.matrix)] <- cancer.colnames[j]
    }
  }

  dev.off()

  write.table(abundance.score.matrix,
    file.path(results_dir, "score_matrix.txt"),
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
  )

  abundance.score.matrix
}
