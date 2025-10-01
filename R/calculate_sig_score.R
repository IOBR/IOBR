#' List of Supported Signature Score Calculation Methods
#'
#' @description
#' A named vector containing the supported methods for calculating signature scores.
#' Methods include PCA, ssGSEA, z-score, and Integration. The names represent display
#' names while values represent internal method identifiers.
#'
#' @format Named character vector with the following elements:
#' \describe{
#'   \item{PCA}{Principal Component Analysis method (internal: "pca")}
#'   \item{ssGSEA}{Single-sample Gene Set Enrichment Analysis (internal: "ssgsea")}
#'   \item{z-score}{Z-score transformation method (internal: "zscore")}
#'   \item{Integration}{Integration method (internal: "integration")}
#' }
#'
#' @export

signature_score_calculation_methods <- c(
  "PCA" = "pca",
  "ssGSEA" = "ssgsea",
  "z-score" = "zscore",
  "Integration" = "integration"
)
##########################################




#' Calculate Signature Score Using PCA Method
#'
#' @description
#' Computes signature scores for gene sets using Principal Component Analysis (PCA).
#' The first principal component is extracted from gene expression data for each
#' signature and used as the signature score.
#'
#' @param pdata Data frame containing phenotype data. If \code{NULL}, a data frame with
#'   \code{Index} and \code{ID} columns is created from \code{eset} column names.
#'   Default is \code{NULL}.
#' @param eset Matrix of normalized gene expression data (CPM, TPM, RPKM, FPKM, etc.)
#'   with genes in rows and samples in columns.
#' @param signature List of gene signatures, where each element is a character vector
#'   of gene names.
#' @param mini_gene_count Integer specifying the minimum number of genes required in
#'   the expression set for a signature to be included. Default is 3.
#' @param column_of_sample Character string specifying the column in \code{pdata}
#'   containing sample identifiers. Default is \code{"ID"}.
#' @param adjust_eset Logical indicating whether to remove features with missing values,
#'   zero standard deviation, or infinite values. Default is \code{FALSE}.
#'
#' @return Tibble containing sample identifiers and signature scores (signatures in
#'   columns, samples in rows). Special combined scores (e.g., TMEscore_CIR) are
#'   calculated if corresponding signature pairs are present.
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD expression data (raw count matrix)
#' data("eset_stad", package = "IOBR")
#' # Transform count data to TPM
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' # Calculate signature scores using PCA method
#' calculate_sig_score_pca(eset = eset, signature = signature_tme)
calculate_sig_score_pca <- function(pdata = NULL,
                                    eset,
                                    signature,
                                    mini_gene_count = 3,
                                    column_of_sample = "ID",
                                    adjust_eset = FALSE) {
  message(paste0("\n", ">>> Calculating signature score using PCA method"))

  # creat pdata if NULL
  if (is.null(pdata)) {
    pdata <- data.frame("Index" = 1:length(colnames(eset)), "ID" = colnames(eset))
  } else {
    pdata <- as.data.frame(pdata)

    if ("ID" %in% colnames(pdata) & !column_of_sample == "ID") {
      colnames(pdata)[which(colnames(pdata) == "ID")] <- "ID2"
      message("In order to prevent duplicate names, the 'ID' column of original pdata was rename into 'ID2' ")
    }

    if (column_of_sample %in% colnames(pdata)) {
      colnames(pdata)[which(colnames(pdata) == column_of_sample)] <- "ID"
    }
  }
  # match phenotype data and gene expression set
  ###########################
  pdata <- pdata[pdata$ID %in% colnames(eset), ]
  eset <- eset[, colnames(eset) %in% pdata$ID]
  eset <- eset[, match(pdata$ID, colnames(eset))]
  ###########################

  # normalization
  if (dim(eset)[2] < 5000) eset <- log2eset(eset = eset)
  if (dim(eset)[2] < 5000) check_eset(eset)
  if (adjust_eset) {
    feas <- feature_manipulation(data = eset, is_matrix = T)
    eset <- eset[rownames(eset) %in% feas, ]
  }

  # eset<-scale(eset,center = T,scale = T)
  ###########################

  # filter signatures
  if (mini_gene_count <= 2) mini_gene_count <- 2
  signature <- signature[lapply(signature, function(x) sum(x %in% rownames(eset) == TRUE)) >= mini_gene_count]
  ###########################

  # calculating signature score
  goi <- names(signature)
  ###########################
  for (sig in goi) {
    pdata[, sig] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop = FALSE]
    pdata[, sig] <- sigScore(tmp, methods = "PCA")
  }
  if ("TMEscoreA_CIR" %in% goi & "TMEscoreB_CIR" %in% goi) {
    pdata[, "TMEscore_CIR"] <- pdata[, "TMEscoreA_CIR"] - pdata[, "TMEscoreB_CIR"]
  }
  if ("TMEscoreA_plus" %in% goi & "TMEscoreB_plus" %in% goi) {
    pdata[, "TMEscore_plus"] <- pdata[, "TMEscoreA_plus"] - pdata[, "TMEscoreB_plus"]
  }

  if ("Index" %in% colnames(pdata)) pdata <- pdata[, -which(colnames(pdata) == "Index")]
  pdata <- tibble::as_tibble(pdata)
  return(pdata)
}
###################################################






#' Calculate Signature Score Using Z-Score Method
#'
#' @description
#' Computes signature scores for gene sets using the z-score transformation method.
#' Gene expression values are standardized and averaged across genes within each
#' signature to produce signature scores.
#'
#' @param pdata Data frame containing phenotype data. If \code{NULL}, a data frame with
#'   \code{Index} and \code{ID} columns is created from \code{eset} column names.
#'   Default is \code{NULL}.
#' @param eset Matrix of normalized gene expression data (CPM, TPM, RPKM, FPKM, etc.)
#'   with genes in rows and samples in columns.
#' @param signature List of gene signatures, where each element is a character vector
#'   of gene names.
#' @param mini_gene_count Integer specifying the minimum number of genes required in
#'   the expression set for a signature to be included. Default is 3.
#' @param column_of_sample Character string specifying the column in \code{pdata}
#'   containing sample identifiers. Default is \code{"ID"}.
#' @param adjust_eset Logical indicating whether to remove features with missing values,
#'   zero standard deviation, or infinite values. Default is \code{FALSE}.
#'
#' @return Data frame containing phenotype data and signature scores (signatures in
#'   columns, samples in rows).
#'
#' @author Dongqiang Zeng
#' @export
#' @examples
#' # Load TCGA-STAD expression data (raw count matrix)
#' data("eset_stad", package = "IOBR")
#' # Transform count data to TPM
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' # Calculate signature scores using z-score method
#' calculate_sig_score_zscore(eset = eset, signature = signature_tme)
calculate_sig_score_zscore <- function(pdata = NULL,
                                       eset,
                                       signature,
                                       mini_gene_count = 3,
                                       column_of_sample = "ID",
                                       adjust_eset = FALSE) {
  message(paste0("\n", ">>> Calculating signature score using z-score method"))

  # creat pdata if NULL
  if (is.null(pdata)) {
    pdata <- data.frame("Index" = 1:length(colnames(eset)), "ID" = colnames(eset))
  } else {
    pdata <- as.data.frame(pdata)

    if ("ID" %in% colnames(pdata) & !column_of_sample == "ID") {
      colnames(pdata)[which(colnames(pdata) == "ID")] <- "ID2"
      message("In order to prevent duplicate names, the 'ID' column of original pdata was rename into 'ID2' ")
    }

    if (column_of_sample %in% colnames(pdata)) {
      colnames(pdata)[which(colnames(pdata) == column_of_sample)] <- "ID"
    }
  }
  # match phenotype data and gene expression set
  ###########################
  pdata <- pdata[pdata$ID %in% colnames(eset), ]
  eset <- eset[, colnames(eset) %in% pdata$ID]
  eset <- eset[, match(pdata$ID, colnames(eset))]
  ###########################
  # normalization
  if (ncol(eset) < 5000) eset <- log2eset(eset = eset)
  if (ncol(eset) < 5000) check_eset(eset)

  if (adjust_eset) {
    feas <- feature_manipulation(data = eset, is_matrix = T)
    eset <- eset[rownames(eset) %in% feas, ]
  }

  # eset<-scale(eset,center = T,scale = T)
  ###########################
  ###########################
  if (mini_gene_count <= 2) mini_gene_count <- 2
  signature <- signature[lapply(signature, function(x) sum(x %in% rownames(eset) == TRUE)) >= mini_gene_count]
  ###########################
  # calculating signature score
  goi <- names(signature)
  for (sig in goi) {
    pdata[, sig] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop = FALSE]
    pdata[, sig] <- sigScore(tmp, methods = "zscore")
  }
  if ("TMEscoreA_CIR" %in% goi & "TMEscoreB_CIR" %in% goi) {
    pdata[, "TMEscore_CIR"] <- pdata[, "TMEscoreA_CIR"] - pdata[, "TMEscoreB_CIR"]
  }
  if ("TMEscoreA_plus" %in% goi & "TMEscoreB_plus" %in% goi) {
    pdata[, "TMEscore_plus"] <- pdata[, "TMEscoreA_plus"] - pdata[, "TMEscoreB_plus"]
  }

  if ("Index" %in% colnames(pdata)) pdata <- pdata[, -which(colnames(pdata) == "Index")]
  pdata <- tibble::as_tibble(pdata)
  return(pdata)
}
###################################################




#' Calculating signature score using ssGSEA method
#'
#' @param pdata phenotype data of input sample;
#' if phenotype data is NULL, create a data frame with `Index` and `ID` contain column names of eset
#' @param eset normalizated  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set; default is 5;
#' the minimal gene count for ssGSEA methods should larger than 5 for the robustness of the calculation
#' @param column_of_sample  Defines in which column of pdata the sample identifier can be found.
#' @param adjust_eset remove variables with missing value, sd =0, and Inf value
#' @param parallel.size default is 1
#'
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @import tibble
#' @author Dongqiang Zeng
#' @examples
#' # Loading TCGA-STAD expresion data(raw count matrix)
#' data("eset_stad", package = "IOBR")
#' # transform count data to tpm
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' # signature score estimation using ssgsea method
#' calculate_sig_score_ssgsea(eset = eset, signature = signature_tme)
#'
calculate_sig_score_ssgsea <- function(pdata = NULL,
                                       eset,
                                       signature,
                                       mini_gene_count = 3,
                                       column_of_sample = "ID",
                                       adjust_eset = FALSE,
                                       parallel.size = 1L) {
  message(paste0("\n", ">>> Calculating signature score using ssGSEA method"))

  signature <- signature[lapply(signature, function(x) sum(x %in% rownames(eset) == TRUE)) >= mini_gene_count]

  if (mini_gene_count <= 5) mini_gene_count <- 5
  ############################

  # creat pdata if NULL
  if (is.null(pdata)) {
    pdata <- data.frame("Index" = 1:length(colnames(eset)), "ID" = colnames(eset))
  } else {
    pdata <- as.data.frame(pdata)

    if ("ID" %in% colnames(pdata) & !column_of_sample == "ID") {
      colnames(pdata)[which(colnames(pdata) == "ID")] <- "ID2"
      message("In order to prevent duplicate names, the 'ID' column of original pdata was rename into 'ID2' ")
    }

    if (column_of_sample %in% colnames(pdata)) {
      colnames(pdata)[which(colnames(pdata) == column_of_sample)] <- "ID"
    }
  }
  # match phenotype data and gene expression set
  ###########################
  pdata <- pdata[pdata$ID %in% colnames(eset), ]
  eset <- eset[, colnames(eset) %in% pdata$ID]

  ##############################
  if (ncol(eset) < 5000) eset <- log2eset(eset = eset)
  if (ncol(eset) < 5000) check_eset(eset)

  if (adjust_eset) {
    feas <- feature_manipulation(data = eset, is_matrix = T)
    eset <- eset[rownames(eset) %in% feas, ]
  }

  ##############################
  # Check the formal argument of GSVA::gsva
  FA <- formals(GSVA::gsva)

  if (is.null(FA[["method"]])) {
    params <- gsvaParam(as.matrix(eset),
      signature,
      minSize = mini_gene_count,
      maxSize = Inf,
      kcdf = "Gaussian",
      tau = 1,
      maxDiff = TRUE,
      absRanking = FALSE
    )

    res <- GSVA::gsva(params,
      verbose = TRUE,
      BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
    )
  } else {
    res <- GSVA::gsva(as.matrix(eset),
      signature,
      method = "ssgsea",
      kcdf = "Gaussian",
      min.sz = mini_gene_count,
      ssgsea.norm = T,
      parallel.sz = parallel.size
    )
  }

  ##############################
  res <- as.data.frame(t(res))
  res <- rownames_to_column(res, var = "ID")

  if ("TMEscoreA_CIR" %in% colnames(res) & "TMEscoreB_CIR" %in% colnames(res)) {
    res[, "TMEscore_CIR"] <- res[, "TMEscoreA_CIR"] - res[, "TMEscoreB_CIR"]
  }

  if ("TMEscoreA_plus" %in% colnames(res) & "TMEscoreB_plus" %in% colnames(res)) {
    res[, "TMEscore_plus"] <- res[, "TMEscoreA_plus"] - res[, "TMEscoreB_plus"]
  }

  pdata <- merge(pdata, res, by = "ID", all.x = F, all.y = F)

  if ("Index" %in% colnames(pdata)) pdata <- pdata[, -which(colnames(pdata) == "Index")]
  pdata <- tibble::as_tibble(pdata)
  return(pdata)
}
###################################################






#' Calculating signature score using Integration method
#'
#' @param pdata phenotype data of input sample;
#' if phenotype data is NULL, create a data frame with `Index` and `ID` contain column names of eset
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set; default is 5;
#' the minimal gene count for ssGSEA methods should larger than 5 for the robustness of the calculation
#' @param adjust_eset default is FALSE, if true, data with Inf or zero sd will be replaced
#' @param parallel.size default is 1
#' @param column_of_sample  Defines in which column of pdata the sample identifier can be found.
#'
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @import tibble
#' @author Dongqiang Zeng
#' @examples
#' # Loading TCGA-STAD expresion data(raw count matrix)
#' data("eset_stad", package = "IOBR")
#' # transform count data to tpm
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' # signature score estimation using PCA, z-score, and ssgsea method
#' calculate_sig_score_integration(eset = eset, signature = signature_tme)
calculate_sig_score_integration <- function(pdata = NULL,
                                            eset,
                                            signature,
                                            mini_gene_count = 2,
                                            column_of_sample = "ID",
                                            adjust_eset = FALSE,
                                            parallel.size = 1L) {
  message(paste0("\n", ">>> Calculating signature score using PCA, z-score and ssGSEA methods"))

  signature <- signature[lapply(signature, function(x) sum(x %in% rownames(eset) == TRUE)) >= mini_gene_count]
  ###########################
  # creat pdata if NULL
  if (is.null(pdata)) {
    pdata <- data.frame("Index" = 1:length(colnames(eset)), "ID" = colnames(eset))
  } else {
    pdata <- as.data.frame(pdata)

    if ("ID" %in% colnames(pdata) & !column_of_sample == "ID") {
      colnames(pdata)[which(colnames(pdata) == "ID")] <- "ID2"
      message("In order to prevent duplicate names, the 'ID' column of original pdata was rename into 'ID2' ")
    }

    if (column_of_sample %in% colnames(pdata)) {
      colnames(pdata)[which(colnames(pdata) == column_of_sample)] <- "ID"
    }
  }
  # match phenotype data and gene expression set
  ###########################
  pdata <- pdata[pdata$ID %in% colnames(eset), ]
  eset <- eset[, colnames(eset) %in% pdata$ID]
  eset <- eset[, match(pdata$ID, colnames(eset))]
  ###########################
  # normalization
  if (ncol(eset) < 5000) eset <- log2eset(eset = eset)
  if (ncol(eset) < 5000) check_eset(eset)

  if (adjust_eset) {
    feas <- feature_manipulation(data = eset, is_matrix = T)
    eset <- eset[rownames(eset) %in% feas, ]
  }
  ##########################
  # eset1<-scale(eset,center = T,scale = T)
  # if(sum(is.na(eset1))>0){
  #   feas<-feature_manipulation(data=eset1,is_matrix = T)
  #   eset1<-eset1[rownames(eset1)%in%feas,]
  # }
  ###########################
  message(paste0("\n", ">>>Step 1: Calculating signature score using PCA method"))
  goi <- names(signature)
  ###########################
  for (sig in goi) {
    pdata[, paste(sig, "_PCA", sep = "")] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop = FALSE]
    pdata[, paste(sig, "_PCA", sep = "")] <- sigScore(tmp, methods = "PCA")
  }
  ############################
  if ("TMEscoreA_CIR" %in% goi & "TMEscoreB_CIR" %in% goi) {
    pdata[, "TMEscore_CIR_PCA"] <- pdata[, "TMEscoreA_CIR_PCA"] - pdata[, "TMEscoreB_CIR_PCA"]
  }
  if ("TMEscoreA_plus" %in% goi & "TMEscoreB_plus" %in% goi) {
    pdata[, "TMEscore_plus_PCA"] <- pdata[, "TMEscoreA_plus_PCA"] - pdata[, "TMEscoreB_plus_PCA"]
  }
  ############################
  message(paste0("\n", ">>>Step 2: Calculating signature score using z-score method"))
  for (sig in goi) {
    pdata[, paste(sig, "_zscore", sep = "")] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop = FALSE]
    pdata[, paste(sig, "_zscore", sep = "")] <- sigScore(tmp, methods = "z-score")
  }
  if ("TMEscoreA_CIR" %in% goi & "TMEscoreB_CIR" %in% goi) {
    pdata[, "TMEscore_CIR_zscore"] <- pdata[, "TMEscoreA_CIR_zscore"] - pdata[, "TMEscoreB_CIR_zscore"]
  }
  if ("TMEscoreA_plus" %in% goi & "TMEscoreB_plus" %in% goi) {
    pdata[, "TMEscore_plus_zscore"] <- pdata[, "TMEscoreA_plus_zscore"] - pdata[, "TMEscoreB_plus_zscore"]
  }
  ############################
  message(paste0("\n", ">>>Step 3: Calculating signature score using ssGSEA method"))

  # Check the formal argument of GSVA::gsva
  FA <- formals(GSVA::gsva)

  if (is.null(FA[["method"]])) {
    params <- gsvaParam(as.matrix(eset),
      signature,
      minSize = mini_gene_count,
      maxSize = Inf,
      kcdf = "Gaussian",
      tau = 1,
      maxDiff = TRUE,
      absRanking = FALSE
    )

    res <- GSVA::gsva(params,
      verbose = TRUE,
      BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
    )
  } else {
    res <- GSVA::gsva(as.matrix(eset),
      signature,
      method = "ssgsea",
      kcdf = "Gaussian",
      min.sz = mini_gene_count,
      ssgsea.norm = T,
      parallel.sz = parallel.size
    )
  }

  #####################
  res <- as.data.frame(t(res))
  if ("TMEscoreA_CIR" %in% colnames(res) & "TMEscoreB_CIR" %in% colnames(res)) {
    res[, "TMEscore_CIR"] <- res[, "TMEscoreA_CIR"] - res[, "TMEscoreB_CIR"]
  }
  if ("TMEscoreA_plus" %in% colnames(res) & "TMEscoreB_plus" %in% colnames(res)) {
    res[, "TMEscore_plus"] <- res[, "TMEscoreA_plus"] - res[, "TMEscoreB_plus"]
  }
  #############################
  colnames(res) <- paste0(colnames(res), "_ssGSEA")
  res <- rownames_to_column(res, var = "ID")

  pdata <- merge(pdata, res, by = "ID", all.x = F, all.y = F)
  pdata <- tibble::as_tibble(pdata)
  return(pdata)
}






#' Calculating signature score on a gene expression dataset
#'
#' @param pdata phenotype data of input sample
#' @param eset normalized  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures;
#' such as `signature_tme`, `signature_metabolism`,`signature_collection`, `go_bp`,`kegg`,`hallmark`
#' @param method he methods currently supported are `pca`, `ssgsea`, `zscore`,`integration`
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set;
#' default is 2 for PCA and z score function
#' @param column_of_sample Defines in which column of the `pdata` the identifier of the sample can be found.
#' @param print_gene_propotion logical, print the proportion of signature genes in gene matrix
#' @param print_filtered_signatures logical, print filtered signatures has gene count less than minimal gene count
#' @param adjust_eset default is FALSE
#' @param parallel.size default is 1
#'
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @author Dongqiang Zeng
#' @examples
#' # Loading TCGA-STAD expresion data(raw count matrix)
#' data("eset_stad", package = "IOBR")
#' # transform count data to tpm
#' eset <- count2tpm(eset_stad, idType = "ensembl")
#' # signature score estimation using PCA, z-score, or ssgsea method
#' calculate_sig_score(eset = eset, signature = signature_tme, method = "pca")
#'
#' @references 1. Hänzelmann S, Castelo R, Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, 7. doi: 10.1186/1471-2105-14-7
#' @references 2. Mariathasan, S., Turley, S., Nickles, D. et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature 554, 544–548 (2018).
#'
calculate_sig_score <- function(pdata = NULL,
                                eset,
                                signature = signature_collection,
                                method = "pca",
                                mini_gene_count = 3,
                                column_of_sample = "ID",
                                print_gene_propotion = FALSE,
                                adjust_eset = FALSE,
                                parallel.size = 1L,
                                print_filtered_signatures = FALSE, ...) {
  if (class(signature) == "list") {
    signature <- lapply(signature, function(x) na.omit(x))
    signature <- lapply(signature, function(x) as.character(x))
    signature <- lapply(signature, function(x) unique(x))
    signature <- lapply(signature, function(x) x[!x == ""])
  }

  ########################################
  if (print_gene_propotion) {
    message(lapply(signature, function(x) summary(x %in% rownames(eset))))
  }
  ########################################
  if (print_filtered_signatures) {
    filtered_signature <- signature[lapply(signature, function(x) sum(x %in% rownames(eset) == TRUE)) <= mini_gene_count]
    message(paste0("\n", " The number of filtered signatures is: ", length(filtered_signature)))
    if (length(filtered_signature) >= 10) {
      message(paste("\n", "10 of signatures with gene count less than ", mini_gene_count, " : "))
      message(paste(" ", names(filtered_signature)[1:10], collapse = "\n"))
      message(paste0(
        "\n", "You can use this function to get all filtered sigantures:", "\n",
        "signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))<= mini_gene_count]"
      ))
    } else if (length(filtered_signature) > 0 & length(filtered_signature) < 10) {
      message(paste(names(filtered_signature)[1:length(filtered_signature)], collapse = "\n"))
    } else if (length(filtered_signature) <= 0) {
      message(paste0("\n", "All signatures are eligible for calculating signature score"))
    }
  }
  ##########################################

  method <- tolower(method)

  if (!method %in% c("zscore", "pca", "ssgsea", "integration")) stop("At present, we only provide three methods to calculate the score: PCA, zscore, ssGSEA")
  # run selected method
  res <- switch(method,
    pca = calculate_sig_score_pca(pdata, eset,
      signature = signature,
      mini_gene_count = mini_gene_count,
      column_of_sample = column_of_sample,
      adjust_eset = adjust_eset, ...
    ),
    ssgsea = calculate_sig_score_ssgsea(pdata, eset,
      signature = signature,
      mini_gene_count = mini_gene_count,
      column_of_sample = column_of_sample,
      adjust_eset = adjust_eset,
      parallel.size = parallel.size, ...
    ),
    zscore = calculate_sig_score_zscore(pdata, eset,
      signature = signature,
      mini_gene_count = mini_gene_count,
      column_of_sample = column_of_sample,
      adjust_eset = adjust_eset, ...
    ),
    integration = calculate_sig_score_integration(pdata, eset,
      signature = signature,
      mini_gene_count = mini_gene_count,
      column_of_sample = column_of_sample,
      adjust_eset = adjust_eset,
      parallel.size = parallel.size, ...
    )
  )
  return(res)
}
