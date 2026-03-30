#' Source code for the TIMER deconvolution method.
#'
#' This code is adapted from https://github.com/hanfeisun/TIMER, which
#' again is an adapted version of the original TIMER source code
#' from http://cistrome.org/TIMER/download.html.
#'
#' The method is described in Li et al. Genome Biology 2016;17(1):174., PMID 27549193.


#' Display Timed Information Messages
#'
#' This function formats and displays informational messages, primarily for timing or logging purposes. It is useful for tracking the progress or stages of execution within scripts, especially in long-running operations. The function outputs messages to the R console with a prefixed '##' to highlight them.
#'
#' @param string A character string representing the message to be displayed.
#'
#' @return None; the function is used for its side effect of printing a message.
#' @export
#' @author Bo Li
#' @examples
#' TimerINFO("Data processing started.")
TimerINFO <- function(string) {
  message(sprintf("## %s\n", string))
}


#' TIMER signatures are cancer specific. This is the list of available cancer types.
#'
#' @export
timer_available_cancers <- c(
  "kich", "blca", "brca", "cesc", "gbm", "hnsc", "kirp", "lgg",
  "lihc", "luad", "lusc", "prad", "sarc", "pcpg", "paad", "tgct",
  "ucec", "ov", "skcm", "dlbc", "kirc", "acc", "meso", "thca",
  "uvm", "ucs", "thym", "esca", "stad", "read", "coad", "chol"
)


#' Remove batch effect of expression set
#'
#' This function removes batch effects between two gene expression data sets, typically representing different sample types such as cancer cells and immune cells. The function aligns the gene names, combines the datasets, applies batch correction using ComBat from the sva package, and then separates the adjusted data sets. Additionally, it aggregates expression data for immune cell types by taking the median expression level across samples.
#'
#' @param cancer.exp A matrix or data frame of cancer cell expression data with genes as rows and samples as columns.
#' @param immune.exp A matrix or data frame of immune cell expression data with genes as rows and samples as columns.
#' @param immune.cellType A vector indicating the cell type for each column in `immune.exp`.
#'
#' @author Bo Li
#'
#' @return A list containing three elements:
#'   - The first element is the batch effect corrected cancer expression data.
#'   - The second element is the batch effect corrected immune expression data.
#'   - The third element is the aggregated immune expression data, with median values calculated for each cell type.
#'
#' @export
#'
#' @examples
#' set.seed(123) # For reproducibility
#' gene_names <- paste0("Gene", 1:100)
#' sample_names_cancer <- paste0("CancerSample", 1:10)
#' cancer.exp <- matrix(runif(1000, 1, 1000),
#'   nrow = 100, ncol = 10,
#'   dimnames = list(gene_names, sample_names_cancer)
#' )
#'
#' # Generate synthetic expression data for immune cells
#' sample_names_immune <- paste0("ImmuneSample", 1:5)
#' immune.exp <- matrix(runif(500, 1, 1000),
#'   nrow = 100, ncol = 5,
#'   dimnames = list(gene_names, sample_names_immune)
#' )
#'
#' # Create a cell type vector for immune samples
#' immune.cellType <- c("T-cell", "B-cell", "T-cell", "NK-cell", "B-cell")
#' names(immune.cellType) <- sample_names_immune
#'
#' result <- RemoveBatchEffect(cancer.exp, immune.exp, immune.cellType)
RemoveBatchEffect <- function(cancer.exp, immune.exp, immune.cellType) {
  ## intersect the gene names of cancer.exp and immune.exp
  tmp.dd <- as.matrix(cancer.exp)
  # tmp <- sapply(strsplit(rownames(cancer.exp), '\\|'),
  #               function(x) x[[1]])
  # rownames(tmp.dd) <- tmp
  # tmp.dd <- as.matrix(tmp.dd[which(nchar(tmp)>1), ])
  tmp.ss <- intersect(rownames(tmp.dd), rownames(immune.exp))

  ## bind cancer and immune expression data into one dataframe
  N1 <- ncol(tmp.dd)

  tmp.dd <- cbind(tmp.dd[tmp.ss, ], immune.exp[tmp.ss, ])
  tmp.dd <- as.matrix(tmp.dd)
  mode(tmp.dd) <- "numeric"

  ## remove batch effects
  N2 <- ncol(immune.exp)
  tmp.batch <- c(rep(1, N1), rep(2, N2))
  rlang::check_installed(c("sva", "BiocParallel"))
  tmp.dd0 <- sva::ComBat(tmp.dd, tmp.batch, c(), BPPARAM = BiocParallel::bpparam("SerialParam"))

  ## separate cancer and immune expression data after batch effect removing
  dd.br <- tmp.dd0[, 1:N1]
  immune.exp.br <- tmp.dd0[, (N1 + 1):(N1 + N2)]

  ## a immune category has multiple samples, use the median expression level for a gene
  tmp0 <- c()
  for (kk in unique(names(immune.cellType))) {
    tmp.vv <- which(names(immune.cellType) == kk)

    if (length(tmp.vv) == 1) {
      median_expression <- apply(immune.exp.br[, tmp.vv, drop = FALSE], 1, median, na.rm = TRUE)
    } else {
      median_expression <- apply(immune.exp.br[, tmp.vv], 1, median, na.rm = TRUE)
    }

    if (!exists("tmp0")) {
      tmp0 <- median_expression
    } else {
      tmp0 <- cbind(tmp0, median_expression)
    }
  }

  immune.exp.agg.br <- tmp0
  colnames(immune.exp.agg.br) <- unique(names(immune.cellType))
  return(list(as.matrix(dd.br), immune.exp.br, immune.exp.agg.br))
}


#' Process Batch Table and Validate Cancer Types
#'
#' This function processes input data containing cancer types and validates
#' each category against a predefined list of supported cancer types
#' (`timer_available_cancers`).
#'
#' The input can be provided in two ways:
#' \itemize{
#'   \item \strong{Batch mode}: A CSV file path specified in \code{args$batch}.
#''         The file is expected to have at least two columns, where the
#'         second column contains cancer categories.
#'   \item \strong{Direct mode}: Vectors provided via \code{args$expression}
#'         and \code{args$category}.
#' }
#'
#' If any cancer category is not recognized, the function stops with an error.
#'
#' @param args A list containing input parameters:
#' \describe{
#'   \item{batch}{Character string. Path to a CSV file (optional).}
#'   \item{expression}{Character vector of expression identifiers (used if \code{batch} is NULL).}
#'   \item{category}{Character vector of cancer types (used if \code{batch} is NULL).}
#' }
#'
#' @return A character matrix with two columns:
#' \itemize{
#'   \item Column 1: expression identifiers
#'   \item Column 2: cancer categories
#' }
#' Each row corresponds to one record.
#'
#' @examples
#' # ---- Batch mode example ----
#' tf <- tempfile(fileext = ".csv")
#' write.table(
#'   data.frame(
#'     "exp1", "lung",
#'     "exp2", "breast"
#'   ),
#'   file = tf,
#'   sep = ",",
#'   row.names = FALSE,
#'   col.names = FALSE,
#'   quote = FALSE
#' )
#'
#' args <- list(batch = tf)
#' # timer_available_cancers <- c("lung", "breast")  # define for testing
#' result <- check_cancer_types(args)
#'
#' # ---- Direct input example ----
#' args <- list(
#'   expression = c("exp1", "exp2"),
#'   category = c("lung", "breast"),
#'   batch = NULL
#' )
#'
#' # timer_available_cancers <- c("lung", "breast")  # define for testing
#' result <- check_cancer_types(args)
check_cancer_types <- function(args) {
  if (!is.null(args$batch)) {
    TimerINFO("Enter batch mode\n")
    cancers <- as.matrix(read.table(args$batch, sep = ",", stringsAsFactors = FALSE))
  } else {
    if (length(args$expression) != length(args$category)) {
      stop("expression and category must have the same length")
    }
    cancers <- cbind(args$expression, args$category)
  }

  for (i in seq_len(nrow(cancers))) {
    cancer.category <- cancers[i, 2]
    if (!(cancer.category %in% timer_available_cancers)) {
      stop(paste("unknown cancers:", cancer.category))
    }
  }

  return(cancers)
}


#' Constrained regression method implemented in Abbas et al., 2009
#'
#' This function implements a constrained regression approach described by Abbas et al. in their 2009 paper.
#' It is designed to estimate the proportions of immune cell types within a mixed cancer tissue sample
#' based on gene expression data. The method iteratively adjusts the regression coefficients to ensure
#' they are non-negative, which corresponds to realistic proportions of cell types.
#'
#' @param XX Matrix representing immune expression data with genes as rows and cell types as columns.
#' @param YY Vector representing cancer expression data with gene expression levels.
#' @param w Optional vector of weights for the regression, default is NA which means no weights are applied.
#'
#' @return A vector with non-negative coefficients representing the proportions of each cell type in the input expression data.
#' @export
#'
#' @examples
#' # Generate some example data
#' XX <- matrix(runif(100), nrow = 10, ncol = 10)
#' colnames(XX) <- paste("CellType", 1:10, sep = "")
#' YY <- runif(10)
#'
#' # Apply the Abbas constrained regression method
#' results <- GetFractions.Abbas(XX, YY)
#' print(results)
GetFractions.Abbas <- function(XX, YY, w = NA) {
  ss.remove <- c()
  ss.names <- colnames(XX)
  while (T) {
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
    if (is.na(w[1])) tmp <- lsfit(tmp.XX, YY, intercept = F) else tmp <- lsfit(tmp.XX, YY, w, intercept = F)
    if (is.null(ncol(tmp.XX))) tmp.beta <- tmp$coefficients[1] else tmp.beta <- tmp$coefficients[1:(ncol(tmp.XX) + 0)]
    if (min(tmp.beta > 0)) break
    ss.remove <- which.min(tmp.beta)
  }
  tmp.F <- rep(0, ncol(XX))
  names(tmp.F) <- colnames(XX)
  tmp.F[ss.names] <- tmp.beta
  return(tmp.F)
}

#' Convert Rowname To Loci
#'
#' This function processes a gene expression data matrix by modifying its row names.
#' It extracts the gene identifier from the row names assuming they contain additional
#' information separated by a '|'. The function is specifically designed to handle
#' row names formatted as 'LOC389332|389332', and it will simplify these to 'LOC389332'.
#'
#' @param cancerGeneExpression A matrix or data frame of cancer gene expression data
#'        with row names in the format 'GENE|ID'. The function will modify these row names
#'        to keep only the gene part before the '|'.
#'
#' @return A matrix containing the modified gene expression data with updated row names.
#'         Rows without a valid identifier (i.e., names not containing '|') are removed.
#'
#' @export
#'
#' @examples
#' # Assume `data` is your original gene expression matrix with compound row names
#' example_data <- matrix(runif(20), ncol = 5)
#' rownames(example_data) <- c("LOC101", "LOC102", "LOC103", "LOC104")
#' # Process the data to convert row names
#' processed_data <- ConvertRownameToLoci(example_data)
#' print(processed_data)
ConvertRownameToLoci <- function(cancerGeneExpression) {
  ## Extract only the loci information for row name

  ## Example of origin row name is 'LOC389332|389332'
  ## Coverted row name is 'LOC389332'

  ## Args:
  ##   geneExpression: the orginal geneExpression load from .Rdata file
  ##
  ## Returns:
  ##   Modified geneExpression

  tmp <- strsplit(rownames(cancerGeneExpression), "\\|")
  tmp <- sapply(tmp, function(x) x[[1]])
  tmp.vv <- which(nchar(tmp) > 1)
  rownames(cancerGeneExpression) <- tmp
  extracted <- as.matrix(cancerGeneExpression[tmp.vv, ])
  colnames(extracted) <- colnames(cancerGeneExpression)
  return(extracted)
}


#' Parse Input Gene Expression Data
#'
#' This function reads gene expression data from a tab-delimited text file,
#' using the first column as row names. The data are converted into a numeric
#' matrix for downstream analysis.
#'
#' @param path Character string. Path to a tab-delimited gene expression file.
#'   The first column should contain gene identifiers, and the remaining columns
#'   should contain expression values for samples.
#'
#' @return A numeric matrix of gene expression values, with genes as rows and
#'   samples as columns.
#'
#' @export
#'
#' @examples
#' # Create a temporary tab-delimited expression file
#' tf <- tempfile(fileext = ".tsv")
#' expr <- data.frame(
#'   gene = c("GeneA", "GeneB", "GeneC"),
#'   Sample1 = c(10, 20, 30),
#'   Sample2 = c(15, 25, 35)
#' )
#' write.table(expr, tf, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' # Parse the gene expression data
#' gene_expression_data <- ParseInputExpression(tf)
#' print(gene_expression_data)
ParseInputExpression <- function(path) {
  ret <- read.delim(path, row.names = 1, check.names = FALSE)
  ret <- as.matrix(ret)
  mode(ret) <- "numeric"
  # ret <- ConvertRownameToLoci(ret)
  return(ret)
}


#' Draw QQ Plot Comparing Cancer and Immune Expression
#'
#' This function creates a quantile-quantile (QQ) plot to compare the gene expression
#' distributions between cancer and immune samples. The QQ plot helps assess if the
#' two distributions are similarly shaped by plotting their quantiles against each other.
#' Points lining up along the diagonal line indicate similar distributions.
#' The function also fits a linear model to the central portion of the data (excluding
#' the lowest 40% and the highest 10% of data points) to highlight the trend.
#'
#' @param cancer.exp A vector of gene expression data for cancer samples.
#' @param immune.exp A vector of gene expression data for immune samples.
#' @param name Optional parameter to add a subtitle with additional information.
#'
#' @return Generates a QQ plot.
#'
#' @export
#'
#' @examples
#' cancer_exp <- rnorm(100, mean = 5, sd = 1.5)
#'
#' immune_exp <- rnorm(100, mean = 5, sd = 1.5)
#'
#' expression_data <- data.frame(
#'   Cancer_Expression = cancer_exp,
#'   Immune_Expression = immune_exp
#' )
#'
#' DrawQQPlot(
#'   cancer.exp = expression_data$Cancer_Expression,
#'   immune.exp = expression_data$Immune_Expression,
#'   name = "Comparison of Gene Expression"
#' )
DrawQQPlot <- function(cancer.exp, immune.exp, name = "") {
  ## q-q plot by sample should look like a straight line.
  ## Extreme values may saturate for Affy array data, but most of the data should align well.
  qq <- qqplot(cancer.exp, immune.exp,
    xlab = "Tumor Expression", ylab = "Ref Expression",
    main = "Sample-Sample Q-Q plot"
  )
  mtext(name, col = "gray11")

  # get part of the points for fit linear, remove bottom 40%, and top 10%
  start <- 0.4 * length(qq$x)
  end <- 0.9 * length(qq$x)
  qq.sub <- list(x = qq$x[start:end], y = qq$y[start:end])
  fit <- lm(y ~ x, data = qq.sub)
  abline(fit, col = "blue")
}

#' Get Outlier Genes
#'
#' This function identifies outlier genes from multiple cancer datasets.
#' It treats the top 5 expressed genes in each sample as outliers and returns a list of unique outlier genes across all samples.
#'
#' @param cancers A dataframe with one column containing paths to gene expression files.
#'
#' @return A vector of unique gene names identified as outliers across all given samples.
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
  ## Return a union of outlier genes.
  ## The top 5 expressed genes in each sample is treated as outlier here.
  outlier.total <- c()
  for (i in seq_len(nrow(cancers))) {
    cancer.expFile <- as.character(cancers[i, 1])
    cancer.expression <- ParseInputExpression(cancer.expFile)
    for (j in seq_len(ncol(cancer.expression))) {
      outlier <- rownames(cancer.expression)[tail(order(cancer.expression[, j]), 5)]
      outlier.total <- c(outlier.total, outlier)
    }
  }
  return(unique(outlier.total))
}


#' Deconvolute Tumor Microenvironment Using TIMER
#'
#' This function performs deconvolution of the tumor microenvironment using the TIMER algorithm. It processes multiple cancer datasets, removes batch effects, and estimates immune cell type abundances. The function relies on specific data and helper functions to manage the expression data and analyze the tumor microenvironment.
#'
#' @param args An environment or list containing parameters and file paths necessary for the function to execute. This should include `outdir` for output directory, `batch` for a file containing paths to expression data and cancer types.
#'
#' @return Returns a matrix of abundance scores for different immune cell types across multiple cancer samples.
#' @export
#'
#' @examples
#' args <- list(
#'   outdir = "path/to/output",
#'   batch = "path/to/batch.csv"
#' )
#' results <- deconvolute_timer.default(args)
deconvolute_timer.default <- function(args) {
  # data("cancer_type_genes")
  # data("immuneCuratedData")

  cancers <- check_cancer_types(args)

  TimerINFO("Loading immune gene expression")

  immune <- immuneCuratedData
  immune.geneExpression <- immune$genes
  immune.cellTypes <- immune$celltypes

  # message(immune$celltypes)

  outlier.genes <- sort(GetOutlierGenes(cancers))


  print(paste("Outlier genes:", paste(outlier.genes, collapse = " ")))

  dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(paste(args$outdir, "/results", sep = ""))) {
    dir.create(paste(args$outdir, "/results", sep = ""))
  }

  abundance.score.matrix <- c()
  pdf(paste(args$outdir, "/results/output.pdf", sep = ""))
  for (i in 1:nrow(cancers)) {
    cancer.expFile <- cancers[i, 1]
    cancer.category <- cancers[i, 2]
    # gene.selected.marker.path <- system.file("extdata", "timer", "precalculated", paste0("genes_", cancer.category, ".RData"),
    #                                          package = "IOBR", mustWork = TRUE)
    cancer.expression <- ParseInputExpression(cancer.expFile)
    index <- !(row.names(cancer.expression) %in% outlier.genes)
    cancer.expression <- cancer.expression[index, , drop = FALSE]
    cancer.colnames <- colnames(cancer.expression)

    TimerINFO(paste("Removing the batch effect of", cancer.expFile))
    for (j in 1:length(cancer.colnames)) {
      DrawQQPlot(cancer.expression[, j], immune.geneExpression[, 1], name = cancer.colnames[j])
    }

    tmp <- RemoveBatchEffect(cancer.expression, immune.geneExpression, immune.cellTypes)
    cancer.expNorm <- tmp[[1]]
    immune.expNormMedian <- tmp[[3]]

    for (j in 1:length(cancer.colnames)) {
      DrawQQPlot(cancer.expNorm[, j], immune.expNormMedian[, 1],
        name = paste("After batch removing and aggregating for", cancer.colnames[j])
      )
    }

    gene.selected.marker <- cancer_type_genes[[which(names(cancer_type_genes) == cancer.category)]]
    gene.selected.marker <- intersect(gene.selected.marker, row.names(cancer.expNorm))
    XX <- immune.expNormMedian[gene.selected.marker, c(-4)]
    YY <- cancer.expNorm[gene.selected.marker, , drop = FALSE]

    for (j in 1:length(cancer.colnames)) {
      fractions <- GetFractions.Abbas(XX, YY[, j])
      # print (paste("Fractions for", cancer.expFile, cancer.colnames[j]))
      # print (fractions)
      barplot(fractions,
        cex.names = 0.8, names.arg = names(fractions), xlab = "cell type", ylab = "abundance",
        main = paste("Abundance estimation for", cancer.colnames[j])
      )
      box()

      abundance.score.matrix <- cbind(abundance.score.matrix, fractions)
      colnames(abundance.score.matrix)[ncol(abundance.score.matrix)] <- cancer.colnames[j]
    }
  }

  dev.off()
  write.table(abundance.score.matrix, paste(args$outdir, "/results/score_matrix.txt", sep = ""),
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
  )

  return(abundance.score.matrix)
}
