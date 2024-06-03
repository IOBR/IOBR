#' Generate Reference Signature Matrix using DESeq2
#'
#' This function performs differential expression analysis using the DESeq2 package
#' to identify genes that are significantly expressed across different cell types,
#' as specified in `pheno`. It calculates median expression levels of these significant
#' genes to form a reference signature matrix.
#'
#' @param dds raw count data from RNA-seq
#' @param pheno  character vector; cell type class of the samples
#' @param FDR numeric; genes with BH adjust p value < FDR are considered significant.
#' @param dat data frame or matrix; normalized transcript quantification data (like FPKM, TPM). Note: cell's median expression level of the identified probes will be the output of reference_matrix.
#' @return A list containing the following elements:
#'         - `reference_matrix`: A dataframe with the median expression values of significant genes across cell types.
#'         - `G`: The optimal number of probes that minimizes the condition number.
#'         - `condition_number`: The minimum condition number corresponding to the optimal number of probes.
#'         - `whole_matrix`: The whole median expression matrix used for further analysis or validation.
#' @export
#' @examples
#' # Example usage assuming 'dds' is a DESeq2 dataset and 'dat' is normalized data:
#' # Load raw count data
#' dds <- matrix(sample(0:1000, 2000, replace = TRUE), nrow = 100, ncol = 20)
#' colnames(dds) <- paste("Sample", 1:20, sep = "_")
#' rownames(dds) <- paste("Gene", 1:100, sep = "_")
#'
#' # Create phenotype data
#' pheno <- rep(c("Type1", "Type2"), each = 10)
#'
#' # Load normalized data
#' dat <- matrix(runif(2000), nrow = 100, ncol = 20)
#' rownames(dat) <- rownames(dds)
#' colnames(dat) <- colnames(dds)
#'
#' # Generate reference signature matrix
#' result <- generateRef_DEseq2(dds = dds, pheno = pheno, FDR = 0.05, dat = dat)
#' print(result$reference_matrix)
generateRef_DEseq2 <- function(dds, pheno, FDR = 0.05, dat){
  if (!all(colnames(dds) == colnames(dat))){
    stop("The sample order of dds and dat much be identical")
  }
  dds <- round(dds, 0)
  keep <- rowSums(dds) >= 0.05 *ncol(dds)
  dds <- dds[keep,]
  dds <- na.omit(dds)
  samples <- data.frame(row.names = colnames(dds), cell = pheno)
  ncells <- unique(pheno)
  res <- list()
  for (i in 1:length(ncells)){
     cellA <- ncells[i]
     samples$group <- NA
     samples$group <- ifelse(samples$cell == cellA, cellA, "others")
     dds2 <- DESeq2::DESeqDataSetFromMatrix(dds,
                                   colData = samples,
                                   design = ~ group)
     dds2 <- DESeq2::DESeq(dds2)
     res[[i]] <- DESeq2::results(dds2, contrast = c("group", cellA, "others"))
    }

  resData <- lapply(res, function(x){
    tmp <- as.data.frame(x[order(x$padj), ])
    tmp <- data.frame(probe = rownames(tmp), tmp)
    tmp <- tmp[tmp$padj < FDR, ]
    return(tmp)
  })
  median_value <- t(dat) %>% as.data.frame() %>%
    split(., pheno) %>% purrr::map(function(x) matrixStats:: colMedians(as.matrix(x)))
  median_value <- do.call(cbind, median_value)
  rownames(median_value) <- rownames(dat)

  con_nums <- c()
  for (i in 50:200){
    probes <- resData %>% purrr::map(function(x)Top_probe(dat = x, i = i)) %>%
      unlist() %>% unique()
    tmpdat <- median_value[probes, ]
    #condition number
    con_num <- kappa(tmpdat)
    con_nums <- c(con_nums, con_num)
  }
  i <- c(50:200)[which.min(con_nums)]
  probes <- resData %>% purrr::map(function(x)Top_probe(dat = x, i = i)) %>%
    unlist() %>% unique()
  reference <- median_value[probes, ]
  reference <- data.frame(NAME = rownames(reference), reference)
  G <- i
  condition_number <- min(con_nums)
  return(list(reference_matrix = reference, G = G, condition_number = condition_number,
              whole_matrix = median_value))
}


#' Top Probe Selector
#'
#' This function extracts the top `i` probes based on their ordering in the provided data frame. 
#' If the number of rows in the data frame is less than or equal to `i`, it returns all probes.
#'
#' @param dat A data frame containing a column named "probe" among other data.
#' @param i An integer indicating the number of top probes to return from the dataset.
#'
#' @return A character vector containing the names of the top `i` probes.
#' @export
#'
#' @examples
#' # Assuming 'dat' is a data frame with at least one column named "probe"
#' dat <- data.frame(probe = c("Probe1", "Probe2", "Probe3", "Probe4", "Probe5"),
#'                   value = c(5, 3, 2, 4, 1))
#' # Get the top 3 probes
#' top_probes <- Top_probe(dat, 3)
#' print(top_probes)
Top_probe <- function(dat, i){
  if (nrow(dat) <= i){
    probe = dat[, "probe"] %>% as.character()
  }else{
    probe = dat[1:i, "probe"] %>% as.character()
  }
  probe
}
