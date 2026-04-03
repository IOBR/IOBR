#' Removing Batch Effect from Expression Sets
#'
#' @description
#' Removes batch effects from expression datasets using sva::ComBat (for
#' microarray/TPM data) or sva::ComBat_seq (for RNA-seq count data). Generates
#' PCA plots to visualize data before and after correction.
#'
#' @param eset1 First expression set (matrix or data frame with genes as rows).
#' @param eset2 Second expression set.
#' @param eset3 Optional third expression set. Use `NULL` if not available.
#' @param id_type Type of gene ID in expression sets (e.g., `"ensembl"`,
#'   `"symbol"`). Required for count data normalization.
#' @param data_type Type of data: `"array"`, `"count"`, or `"tpm"`.
#'   Default is `"array"`.
#' @param cols Color scale for PCA plot. Default is `"normal"`.
#' @param palette Color palette for PCA plot. Default is `"jama"`.
#' @param log2 Whether to perform log2 transformation. Default is `TRUE`.
#'   Ignored for count data.
#' @param check_eset Whether to check expression sets for errors.
#'   Default is `TRUE`.
#' @param adjust_eset Whether to adjust expression sets by removing problematic
#'   features. Default is `TRUE`.
#' @param repel Whether to add repelling labels to PCA plot. Default is `FALSE`.
#' @param path Directory where results should be saved. Default is `NULL`
#'   (display only).
#'
#' @return Expression matrix after batch correction.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @references
#' Zhang Y, et al. ComBat-seq: batch effect adjustment for RNA-seq count data.
#' NAR Genomics and Bioinformatics. 2020;2(3):lqaa078.
#' doi:10.1093/nargab/lqaa078
#'
#' Leek JT, et al. The sva package for removing batch effects and other
#' unwanted variation in high-throughput experiments. Bioinformatics.
#' 2012;28(6):882-883.
#'
#' @examples
#' eset_stad <- load_data("eset_stad")
#' eset_blca <- load_data("eset_blca")
#' eset_corrected <- remove_batcheffect(
#'   eset_stad[1:1000, 1:5], eset_blca[1:1000, 1:5],
#'   id_type = "ensembl",
#'   data_type = "count"
#' )
remove_batcheffect <- function(eset1,
                               eset2,
                               eset3 = NULL,
                               id_type,
                               data_type = c("array", "count", "tpm"),
                               cols = "normal",
                               palette = "jama",
                               log2 = TRUE,
                               check_eset = TRUE,
                               adjust_eset = TRUE,
                               repel = FALSE,
                               path = NULL) {
  # TODO
  # Ignoring unknown labels:
  # linetype : "batch"

  # Validate data_type
  data_type <- rlang::arg_match(data_type)

  # Create output folder if path provided
  if (!is.null(path)) {
    path <- creat_folder(path)
  }

  # Log2 transformation for non-count data
  if (data_type != "count" && log2) {
    eset1 <- log2eset(eset1)
    eset2 <- log2eset(eset2)
    if (!is.null(eset3)) {
      eset3 <- log2eset(eset3)
    }
  }

  # Check expression sets
  if (check_eset) {
    check_eset(eset1)
    check_eset(eset2)
    if (!is.null(eset3)) check_eset(eset3)
  }

  # Adjust expression sets
  if (adjust_eset) {
    eset1 <- eset1[rownames(eset1) %in%
      feature_manipulation(eset1, is_matrix = TRUE), ]
    eset2 <- eset2[rownames(eset2) %in%
      feature_manipulation(eset2, is_matrix = TRUE), ]
    if (!is.null(eset3)) {
      eset3 <- eset3[rownames(eset3) %in%
        feature_manipulation(eset3, is_matrix = TRUE), ]
    }
  }

  # Find common genes
  if (is.null(eset3)) {
    comgene <- intersect(rownames(eset1), rownames(eset2))
    comgene <- comgene[nzchar(comgene) & !is.na(comgene)]

    cli::cli_alert_info(
      "The two expression matrices share {length(comgene)} features"
    )

    combined.expr <- cbind(
      eset1[comgene, , drop = FALSE],
      eset2[comgene, , drop = FALSE]
    )
    batch <- data.frame(
      ID = colnames(combined.expr),
      batch = rep(c("eset1", "eset2"), times = c(ncol(eset1), ncol(eset2)))
    )
  } else {
    comgene <- intersect(
      intersect(rownames(eset1), rownames(eset2)),
      rownames(eset3)
    )
    comgene <- comgene[nzchar(comgene) & !is.na(comgene)]

    cli::cli_alert_info(
      "The three expression matrices share {length(comgene)} features"
    )

    combined.expr <- cbind(
      eset1[comgene, , drop = FALSE],
      eset2[comgene, , drop = FALSE],
      eset3[comgene, , drop = FALSE]
    )
    batch <- data.frame(
      ID = colnames(combined.expr),
      batch = rep(c("eset1", "eset2", "eset3"),
        times = c(ncol(eset1), ncol(eset2), ncol(eset3))
      )
    )
  }

  # Batch correction
  rlang::check_installed("sva")

  if (data_type != "count") {
    cli::cli_alert_info("Processing method: sva::ComBat")
    modcombat <- stats::model.matrix(~1, data = batch)
    combined.expr.combat <- sva::ComBat(
      dat = as.matrix(combined.expr),
      batch = batch$batch,
      mod = modcombat
    )
    rlang::check_installed("preprocessCore")
    combined.expr.combat <- preprocessCore::normalize.quantiles(
      as.matrix(combined.expr.combat),
      keep.names = TRUE
    )
    prefix <- data_type
  } else {
    cli::cli_alert_info("Processing method: sva::ComBat_seq")
    combined.expr.combat <- sva::ComBat_seq(
      as.matrix(combined.expr),
      batch = batch$batch
    )

    # Convert corrected counts to TPM
    eset2_tpm <- count2tpm(
      countMat = combined.expr.combat,
      idType = id_type,
      source = "local"
    )
    eset2_tpm <- log2eset(eset2_tpm)
    cli::cli_alert_info("Count data processed with ComBat_seq")
    prefix <- "count"
  }

  # Generate PCA plots if FactoMineR is available
  if (requireNamespace("FactoMineR", quietly = TRUE) &&
      requireNamespace("factoextra", quietly = TRUE)) {
    p1 <- iobr_pca(
      data = combined.expr,
      is.matrix = TRUE,
      scale = TRUE,
      is.log = TRUE,
      pdata = batch,
      id_pdata = "ID",
      group = "batch",
      cols = cols,
      palette = palette,
      repel = repel,
      ncp = 3,
      axes = c(1, 2),
      addEllipses = TRUE
    ) + ggplot2::ggtitle(paste("Before correction:", prefix))

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
      ncp = 3,
      axes = c(1, 2),
      addEllipses = TRUE
    ) + ggplot2::ggtitle(paste("After correction:", prefix))

    if (data_type == "count") {
      p3 <- iobr_pca(
        data = eset2_tpm,
        is.matrix = TRUE,
        scale = TRUE,
        is.log = TRUE,
        pdata = batch,
        id_pdata = "ID",
        group = "batch",
        cols = cols,
        palette = palette,
        repel = repel,
        ncp = 3,
        axes = c(1, 2),
        addEllipses = TRUE
      ) + ggplot2::ggtitle("After correction: count2TPM")

      p <- patchwork::wrap_plots(p1, p2, p3, nrow = 1)
    } else {
      p <- patchwork::wrap_plots(p1, p2, nrow = 1)
    }

    print(p)

    # Save plot if path provided
    if (!is.null(path)) {
      num_plots <- if (data_type == "count") 3 else 2
      width <- num_plots * 5

      ggplot2::ggsave(
        p,
        filename = paste0("0-PCA-of-", num_plots, "-eset.pdf"),
        width = width,
        height = 5,
        path = path$folder_name
      )
    }
  }

  invisible(combined.expr.combat)
}
