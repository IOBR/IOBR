#' Calculate Immunophenoscore (IPS)
#'
#' @description
#' Calculates Immunophenoscore (IPS) from gene expression data. IPS is a
#' composite score measuring immunophenotype based on four major categories:
#' MHC molecules, immunomodulators, effector cells, and suppressor cells.
#'
#' @param eset Gene expression matrix with official human gene symbols (HGNC)
#'   as rownames. Expression values should be log2(TPM+1) or will be transformed
#'   if max value > 100.
#' @param plot Logical. Whether to generate immunophenogram plots. Default is FALSE.
#' @param project Character string for project identifier. Default is NULL.
#'
#' @return Data frame containing:
#' \describe{
#'   \item{MHC}{MHC molecules score}
#'   \item{EC}{Effector cells score}
#'   \item{SC}{Suppressor cells score}
#'   \item{CP}{Checkpoints/Immunomodulators score}
#'   \item{AZ}{Aggregate score (sum of MHC, CP, EC, SC)}
#'   \item{IPS}{Immunophenoscore (0-10 scale)}
#' }
#'
#' @export
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import purrr
#' @import stringr
#' @import ggplot2
#' @import grid
#' @examples
#' if (interactive()) {
#'   example_genes <- c(
#'     "HLA-A", "HLA-B", "HLA-C", "CD274", "PDCD1", "CTLA4",
#'     "CD8A", "CD8B", "GZMB", "PRF1", "FOXP3", "IL10"
#'   )
#'   sim_eset <- as.data.frame(matrix(
#'     rnorm(length(example_genes) * 5, mean = 5, sd = 2),
#'     nrow = length(example_genes), ncol = 5
#'   ))
#'   rownames(sim_eset) <- example_genes
#'   colnames(sim_eset) <- paste0("Sample", 1:5)
#'   ips_result <- IPS_calculation(eset = sim_eset, project = "Example", plot = FALSE)
#'   if (!is.null(ips_result)) head(ips_result)
#' }
IPS_calculation <- function(eset, project = NULL, plot = FALSE) {
  if (is.null(eset)) return(NULL)

  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort("{.arg eset} must be a matrix or data frame")
  }
  if (ncol(eset) == 0) {
    cli::cli_abort("{.arg eset} must have at least one sample")
  }

  if (max(eset, na.rm = TRUE) > 100) {
    eset <- log2(eset + 1)
  }

  gene_expression <- eset
  sample_names <- colnames(gene_expression)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample", seq_len(ncol(gene_expression)))
    colnames(gene_expression) <- sample_names
  }

  IPSG <- load_data("ips_gene_set")
  if (is.null(IPSG)) return(NULL)

  IPSG <- IPSG[IPSG$GENE %in% rownames(gene_expression), ]

  if (nrow(IPSG) == 0) {
    cli::cli_abort("No IPS genes found in the expression data")
  }

  unique_ips_genes <- unique(IPSG$NAME)

  GVEC <- row.names(gene_expression)
  VEC <- as.vector(IPSG$GENE)
  ind <- which(is.na(match(VEC, GVEC)))
  MISSING_GENES <- VEC[ind]

  if (length(MISSING_GENES) > 0) {
    cli::cli_alert_warning("Differently named or missing genes: {paste(MISSING_GENES, collapse = ', ')}")
  }

  n_samples <- ncol(eset)
  MHC <- numeric(n_samples)
  CP <- numeric(n_samples)
  EC <- numeric(n_samples)
  SC <- numeric(n_samples)
  AZ <- numeric(n_samples)
  IPS <- numeric(n_samples)

  for (i in seq_len(n_samples)) {
    GE <- gene_expression[, i]
    mGE <- mean(GE, na.rm = TRUE)
    sGE <- sd(GE, na.rm = TRUE)
    Z1 <- (gene_expression[as.vector(IPSG$GENE), i] - mGE) / sGE
    W1 <- IPSG$WEIGHT

    MIG <- numeric(length(unique_ips_genes))
    WEIGHT <- numeric(length(unique_ips_genes))

    for (k in seq_along(unique_ips_genes)) {
      gen <- unique_ips_genes[k]
      MIG[k] <- mean(Z1[which(as.vector(IPSG$NAME) == gen)], na.rm = TRUE)
      WEIGHT[k] <- mean(W1[which(as.vector(IPSG$NAME) == gen)], na.rm = TRUE)
    }

    WG <- MIG * WEIGHT
    MHC[i] <- mean(WG[seq_len(10)], na.rm = TRUE)
    CP[i] <- mean(WG[11:20], na.rm = TRUE)
    EC[i] <- mean(WG[21:24], na.rm = TRUE)
    SC[i] <- mean(WG[25:26], na.rm = TRUE)
    AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i], na.rm = TRUE)
    IPS[i] <- ipsmap(AZ[i])

    if (plot) {
      file_name <- file.path(tempdir(), "IPS-Results")
      if (!dir.exists(file_name)) dir.create(file_name, recursive = TRUE)

      plotpath <- file.path(file_name, "IPS_plot_results")
      if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)

      data_a <- data.frame(
        start = c(0, 2.5, 5, 7.5, 10, 15, seq(20, 39), 0, 10, 20, 30),
        end = c(2.5, 5, 7.5, 10, 15, seq(20, 40), 10, 20, 30, 40),
        y1 = c(rep(2.6, 26), rep(0.4, 4)),
        y2 = c(rep(5.6, 26), rep(2.2, 4)),
        z = c(MIG[c(21:26, 11:20, 1:10)], EC[i], SC[i], CP[i], MHC[i]),
        vcol = c(
          unlist(lapply(MIG[c(21:26, 11:20, 1:10)], mapcolors)),
          unlist(lapply(c(EC[i], SC[i], CP[i], MHC[i]), mapbw))
        ),
        label = c(unique_ips_genes[c(21:26, 11:20, 1:10)], "EC", "SC", "CP", "MHC")
      )
      data_a$label <- factor(data_a$label, levels = unique(data_a$label))

      plot_a1 <- ggplot2::ggplot() +
        ggplot2::geom_rect(
          data = data_a, mapping = ggplot2::aes(xmin = start, xmax = end, ymin = y1, ymax = y2, fill = label),
          size = 0.5, color = "black", alpha = 1
        ) +
        ggplot2::coord_polar() +
        ggplot2::scale_y_continuous(limits = c(0, 6)) +
        ggplot2::scale_fill_manual(values = as.vector(data_a$vcol), guide = FALSE) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.margin = grid::unit(0, "mm"), panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "white"),
          axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()
        ) +
        ggplot2::geom_text(ggplot2::aes(x = 5, y = 1.3, label = "EC"), size = 4) +
        ggplot2::geom_text(ggplot2::aes(x = 15, y = 1.3, label = "SC"), size = 4) +
        ggplot2::geom_text(ggplot2::aes(x = 25, y = 1.3, label = "CP"), size = 4) +
        ggplot2::geom_text(ggplot2::aes(x = 35, y = 1.3, label = "MHC"), size = 4)

      plot_a2 <- plot_a1 +
        geom_text(aes(x = 1.25, y = 4.1, label = "+ Act CD4"), angle = 78.75, size = 4) +
        geom_text(aes(x = 3.75, y = 4.1, label = "+ Act CD8"), angle = 56.25, size = 4) +
        geom_text(aes(x = 6.25, y = 4.1, label = "+ Tem CD4"), angle = 33.75, size = 4) +
        geom_text(aes(x = 8.75, y = 4.1, label = "+ Tem CD8"), angle = 11.25, size = 4) +
        geom_text(aes(x = 17.5, y = 4.1, label = "- MDSC"), angle = -67.5, size = 4) +
        geom_text(aes(x = 12.5, y = 4.1, label = "- Treg"), angle = -22.5, size = 4)

      plot_a3 <- plot_a2 +
        geom_text(aes(x = 20.5, y = 4.1, label = "PD-1 -"), angle = 85.5, size = 4) +
        geom_text(aes(x = 21.5, y = 4.1, label = "CTLA4 -"), angle = 76.5, size = 4) +
        geom_text(aes(x = 22.5, y = 4.1, label = "LAG3 -"), angle = 67.5, size = 4) +
        geom_text(aes(x = 23.5, y = 4.1, label = "TIGIT -"), angle = 58.5, size = 4) +
        geom_text(aes(x = 24.5, y = 4.1, label = "TIM3 -"), angle = 49.5, size = 4) +
        geom_text(aes(x = 25.5, y = 4.1, label = "PD-L1 -"), angle = 40.5, size = 4) +
        geom_text(aes(x = 26.5, y = 4.1, label = "PD-L2 -"), angle = 31.5, size = 4) +
        geom_text(aes(x = 27.5, y = 4.1, label = "CD27 +"), angle = 22.5, size = 4) +
        geom_text(aes(x = 28.5, y = 4.1, label = "ICOS +"), angle = 13.5, size = 4) +
        geom_text(aes(x = 29.5, y = 4.1, label = "IDO1 -"), angle = 4.5, size = 4)

      plot_a4 <- plot_a3 +
        geom_text(aes(x = 30.5, y = 4.1, label = "B2M +"), angle = -4.5, size = 4) +
        geom_text(aes(x = 31.5, y = 4.1, label = "TAP1 +"), angle = -13.5, size = 4) +
        geom_text(aes(x = 32.5, y = 4.1, label = "TAP2 +"), angle = -22.5, size = 4) +
        geom_text(aes(x = 33.5, y = 4.1, label = "HLA-A +"), angle = -31.5, size = 4) +
        geom_text(aes(x = 34.5, y = 4.1, label = "HLA-B +"), angle = -40.5, size = 4) +
        geom_text(aes(x = 35.5, y = 4.1, label = "HLA-C +"), angle = -49.5, size = 4) +
        geom_text(aes(x = 36.5, y = 4.1, label = "HLA-DPA1 +"), angle = -58.5, size = 4) +
        geom_text(aes(x = 37.5, y = 4.1, label = "HLA-DPB1 +"), angle = -67.5, size = 4) +
        geom_text(aes(x = 38.5, y = 4.1, label = "HLA-E +"), angle = -76.5, size = 4) +
        geom_text(aes(x = 39.5, y = 4.1, label = "HLA-F +"), angle = -85.5, size = 4)

      plot_a5 <- plot_a4 +
        ggplot2::geom_text(ggplot2::aes(x = 0, y = 6, label = paste("Immunophenoscore: ", IPS[i], sep = "")), angle = 0, size = 6, vjust = -0.5) +
        ggplot2::theme(axis.title = ggplot2::element_blank())

      plot_a <- plot_a5 +
        ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "mm")) +
        ggplot2::geom_text(vjust = 1.15, hjust = 0, ggplot2::aes(x = 25.5, y = 6, label = "\n\n\n\n   MHC: Antigen Processing                                 EC: Effector Cells\n   CP: Checkpoints | Immunomodulators              SC: Suppressor Cells\n\n", hjust = 0), size = 4)

      data_b <- data.frame(
        start = rep(0, 23), end = rep(0.7, 23), y1 = seq(0, 22, by = 1), y2 = seq(1, 23, by = 1),
        z = seq(-3, 3, by = 6 / 22), vcol = c(unlist(lapply(seq(-3, 3, by = 6 / 22), mapcolors))),
        label = LETTERS[1:23]
      )
      data_b_ticks <- data.frame(x = rep(1.2, 7), value = seq(-3, 3, by = 1), y = seq(0, 6, by = 1) * (22 / 6) + 0.5)

      legendtheme <- ggplot2::theme(
        plot.margin = grid::unit(c(2, 0, 2, 0), "inch"), panel.margin = grid::unit(0, "null"),
        panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "white"), axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank()
      )

      plot_b <- ggplot2::ggplot(hjust = 0) +
        ggplot2::geom_rect(data = data_b, mapping = ggplot2::aes(xmin = start, xmax = end, ymin = y1, ymax = y2, fill = label), size = 0.5, color = "black", alpha = 1) +
        ggplot2::scale_x_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
        ggplot2::scale_fill_manual(values = as.vector(data_b$vcol), guide = FALSE) +
        ggplot2::geom_text(data = data_b_ticks, ggplot2::aes(x = x, y = y, label = value), hjust = "inward", size = 4) +
        ggplot2::theme_bw() +
        legendtheme +
        ggplot2::ylab("Sample-wise (averaged) z-score")

      data_c <- data.frame(
        start = rep(0, 23), end = rep(0.7, 23), y1 = seq(0, 22, by = 1), y2 = seq(1, 23, by = 1),
        z = seq(-2, 2, by = 4 / 22), vcol = c(unlist(lapply(seq(-2, 2, by = 4 / 22), mapbw))),
        label = LETTERS[1:23]
      )
      data_c_ticks <- data.frame(x = rep(1.2, 5), value = seq(-2, 2, by = 1), y = seq(0, 4, by = 1) * (22 / 4) + 0.5)

      plot_c <- ggplot2::ggplot() +
        ggplot2::geom_rect(
          data = data_c, mapping = ggplot2::aes(xmin = start, xmax = end, ymin = y1, ymax = y2, fill = label),
          size = 0.5, color = "black", alpha = 1
        ) +
        ggplot2::scale_x_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
        ggplot2::scale_fill_manual(values = as.vector(data_c$vcol), guide = FALSE) +
        ggplot2::geom_text(data = data_c_ticks, ggplot2::aes(x = x, y = y, label = value), hjust = "inward", size = 4) +
        ggplot2::theme_bw() +
        legendtheme +
        ggplot2::ylab("Weighted z-score")

      filename <- file.path(plotpath, paste0("IPS_", sample_names[i], ".pdf"))
      pdf(filename, width = 10, height = 8)
      gridExtra::grid.arrange(plot_a, plot_b, plot_c, ncol = 3, widths = c(0.8, 0.1, 0.1))
      dev.off()
    }
  }

  res <- data.frame(ID = sample_names, MHC = MHC, EC = EC, SC = SC, CP = CP, AZ = AZ, IPS = IPS)

  if (!is.null(project)) {
    res$ProjectID <- project
    res <- res[, c(ncol(res), seq_len(ncol(res) - 1))]
  }

  res <- tibble::column_to_rownames(res, var = "ID")

  return(res)
}
######################################
######################################
######################################
