#' Removing batch effect of two or three expression set
#'
#' @description
#' This function is designed to remove batch effects from given expression datasets and visualize the corrected data using principal component analysis (PCA). It takes three expression datasets as input and performs batch effect correction using the "sva::ComBat" or "sva::ComBat_seq" methods. The function then generates PCA plots to compare the data before and after correction. The PCA plots are customized based on the specified parameters like data type, color palette, log transformation, and path for saving the plots.
#'
#' @param eset1 These are the expression sets for which you want to remove the batch effect.
#' @param eset2 These are the expression sets for which you want to remove the batch effect.
#' @param eset3 These are the expression sets for which you want to remove the batch effect. Input 'NULL' for eset3 if not available.
#' @param check_eset Logical, determining whether to check the eset or not. If TRUE, checks for errors in the expression set. Default is TRUE.
#' @param palette Color palette to use for the plot. Default is 'jama'.
#' @param log2  it performs log2 transformation of the data. Defaults to TRUE.
#' @param path Directory where the results should be saved. If it is NULL, results will be saved in a folder "Combat_PCA".
#' @param adjust_eset Logical, whether to adjust the expression set by manipulating the features. Default is TRUE.
#' @param id_type Type of id present in the expression set (like Ensembl ID or Gene Symbol).
#' @param data_type Type of data in the expression set ("array", "count", or "tpm"). Default is "array".
#' @param cols Color scale to use for the PCA plot. Default is 'normal'.
#' @param repel whether to add labels to the PCA plot. Default is FALSE.
#'
#' @return eset after batch correction
#' @export
#' @author Dongqiang Zeng
#' @references Yuqing Zhang and others, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078
#' @references Leek, J. T., Johnson, W. E., Parker, H. S., Jaffe, A. E., & Storey, J. D. (2012). The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics, 28(6), 882-883.
#' @examples
#' data("eset_stad", package = "IOBR")
#' data("eset_blca", package = "IOBR")
#' # The returned matrix is the count matrix after removing the batches.
#' eset <- remove_batcheffect(eset_stad, eset_blca, id_type = "ensembl", data_type = "count")
remove_batcheffect <- function(eset1, eset2, eset3 = NULL, id_type, data_type = c("array", "count", "tpm"), cols = "normal", palette = "jama",
                               log2 = TRUE, check_eset = TRUE, adjust_eset = TRUE, repel = FALSE, path = "Result_PCA") {
  if (is.null(path)) path <- "Combat_PCA"
  path <- creat_folder(path)
  #######################################################
  if (!data_type == "count") {
    if (log2) {
      eset2 <- log2eset(eset2)
      eset1 <- log2eset(eset1)
    }
  }
  ###########################
  if (check_eset) {
    check_eset(eset1)
    check_eset(eset2)
  }
  ###########################
  if (adjust_eset) {
    feas <- feature_manipulation(data = eset1, is_matrix = T)
    eset1 <- eset1[rownames(eset1) %in% feas, ]

    feas <- feature_manipulation(data = eset2, is_matrix = T)
    eset2 <- eset2[rownames(eset2) %in% feas, ]
  }
  ###########################

  if (!is.null(eset3)) {
    if (log2) eset3 <- log2eset(eset3)
    if (check_eset) check_eset(eset3)
    if (adjust_eset) {
      feas <- feature_manipulation(data = eset3, is_matrix = T)
      eset3 <- eset3[rownames(eset3) %in% feas, ]
    }
  }

  if (is.null(eset3)) {
    comgene <- intersect(rownames(eset1), rownames(eset2))
    comgene <- comgene[!comgene == ""]
    comgene <- comgene[!is.na(comgene)]

    message(paste0(">>>== The two expression matrices share ", length(comgene), " features in common. "))
    combined.expr <- cbind.data.frame(
      eset1[comgene, ],
      eset2[comgene, ]
    )
    batch <- data.frame("ID" = colnames(combined.expr), "batch" = rep(c("eset1", "eset2"), times = c(ncol(eset1), ncol(eset2))))
  }
  ######################

  if (!is.null(eset3)) {
    comgene <- intersect(intersect(rownames(eset1), rownames(eset2)), rownames(eset3))
    comgene <- comgene[!comgene == ""]
    comgene <- comgene[!is.na(comgene)]

    message(paste0(">>>== The three expression matrices share ", length(comgene), " features in common. "))
    combined.expr <- cbind.data.frame(
      eset1[comgene, ],
      eset2[comgene, ],
      eset3[comgene, ]
    )
    batch <- data.frame("ID" = colnames(combined.expr), "batch" = rep(c("eset1", "eset2", "eset3"), times = c(ncol(eset1), ncol(eset2), ncol(eset3))))
  }
  ########################################################################

  if (!data_type == "count") {
    message(">>>=== Processing method: sva:: ComBat")
    modcombat <- model.matrix(~1, data = batch)
    combined.expr.combat <- as.data.frame(sva::ComBat(dat = as.matrix(combined.expr), batch = batch$batch, mod = modcombat))
    combined.expr.combat <- preprocessCore::normalize.quantiles(as.matrix(combined.expr.combat), keep.names = TRUE)
    prefix <- data_type
  } else if (data_type == "count") {
    message(">>>=== Processing method: sva:: ComBat_seq")
    combined.expr.combat <- sva::ComBat_seq(as.matrix(combined.expr), batch = batch$batch)
    # print(head(combined.expr.combat))
    eset2_tpm <- count2tpm(countMat = combined.expr.combat, idType = id_type, source = "local")
    eset2_tpm <- log2eset(eset2_tpm)
    message(">>>=== Count data after proccessing sva::ComBat_seq will be return...")
    prefix <- "count"
  }

  ########################################################################
  p1 <- iobr_pca(
    data = combined.expr,
    is.matrix = TRUE,
    scale = TRUE,
    is.log = TRUE,
    # geom.ind  = "point",
    pdata = batch,
    id_pdata = "ID",
    group = "batch",
    cols = cols,
    palette = palette,
    repel = repel,
    ncp = 3,
    axes = c(1, 2),
    addEllipses = TRUE
  )

  p2 <- iobr_pca(
    data = combined.expr.combat,
    is.matrix = TRUE,
    scale = TRUE,
    is.log = TRUE,
    pdata = batch,
    id_pdata = "ID",
    group = "batch",
    cols = cols,
    palette = palette,
    repel = repel,
    ncp = 3, axes = c(1, 2),
    addEllipses = TRUE
  )
  p1 <- p1 + ggtitle(paste0("Data before correction: ", prefix))
  p2 <- p2 + ggtitle(paste0("Data after correction: ", prefix))

  if (data_type == "count") {
    p3 <- iobr_pca(
      data = eset2_tpm,
      is.matrix = TRUE,
      scale = TRUE, is.log = TRUE, # 取对数+scale
      pdata = batch, id_pdata = "ID", group = "batch",
      cols = cols,
      palette = palette,
      repel = repel,
      ncp = 3, axes = c(1, 2),
      addEllipses = TRUE
    )
    p3 <- p3 + ggtitle("Data after correction: count2TPM")

    p <- p1 | p2 | p3
    # combined.expr.combat <- eset2_tpm
  } else {
    p <- p1 | p2
  }

  print(p)
  ########################################
  if (!data_type == "count") {
    num <- 2
    width <- 2 * 5
  } else {
    num <- 3
    width <- 3 * 5
  }
  ggsave(p, filename = paste0("0-PCA-of-", num, "-eset.pdf"), width = width, height = 5, path = path$folder_name)
  ########################################

  return(combined.expr.combat)
}
