#' Identify Outlier Samples in Gene Expression Data
#'
#' @description
#' Analyzes gene expression data to identify potential outlier samples using
#' connectivity analysis via the WGCNA package. Calculates normalized
#' adjacency and connectivity z-scores for each sample, generates connectivity
#' plots, and optionally performs hierarchical clustering.
#'
#' @param eset Numeric matrix. Gene expression data with genes as rows and
#'   samples as columns.
#' @param yinter Numeric. Z-score threshold for identifying outliers.
#'   Default is -3.
#' @param project Character or `NULL`. Output directory path for saving plots.
#'   Required if `save = TRUE`. Default is `NULL`.
#' @param plot_hculst Logical. Whether to plot hierarchical clustering.
#'   Default is `FALSE`.
#' @param show_plot Logical. Whether to display the connectivity plot.
#'   Default is `TRUE`.
#' @param index Integer or `NULL`. Index for output file naming.
#'   Default is `NULL`.
#' @param save Logical. Whether to save plots to files. Default is `FALSE`.
#'
#' @return Character vector of sample names identified as potential outliers.
#'
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' \donttest{
#' eset_tme_stad <- load_data("eset_tme_stad")
#' outs <- find_outlier_samples(eset = eset_tme_stad)
#' print(outs)
#' }
find_outlier_samples <- function(eset, yinter = -3, project = NULL,
                                 plot_hculst = FALSE, show_plot = TRUE,
                                 index = NULL, save = FALSE) {
  rlang::check_installed("WGCNA")

  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort("{.arg eset} must be a matrix or data frame")
  }
  if (nrow(eset) < 2 || ncol(eset) < 2) {
    cli::cli_abort("{.arg eset} must have at least 2 rows and 2 columns")
  }
  if (!is.numeric(yinter) || length(yinter) != 1) {
    cli::cli_abort("{.arg yinter} must be a single numeric value")
  }

  if (save && is.null(project)) {
    cli::cli_abort("{.arg project} must be provided when save = TRUE")
  }

  path <- if (save) creat_folder(project) else NULL
  if (is.null(index)) index <- 1

  if (save && plot_hculst) {
    tree.combat <- eset %>%
      t() %>%
      stats::dist() %>%
      stats::hclust(method = "average")

    grDevices::pdf(
      paste0(path$abspath, index, "-1-clusteringplot.pdf"),
      width = 20, height = 10
    )
    plot(tree.combat, main = paste0("1-", "Hierarchical Clustering Samples"))
    grDevices::dev.off()
  }

  normalized.adjacency <- (0.5 + 0.5 * WGCNA::bicor(eset))^2
  network.summary <- WGCNA::fundamentalNetworkConcepts(normalized.adjacency)
  connectivity <- network.summary$Connectivity
  connectivity.zscore <- (connectivity - mean(connectivity)) / sqrt(stats::var(connectivity))

  connectivity.plot <- data.frame(
    Sample.Name = names(connectivity.zscore),
    Z.score = connectivity.zscore,
    Sample.Num = seq_along(connectivity.zscore)
  )

  p <- ggplot2::ggplot(
    connectivity.plot,
    ggplot2::aes(x = .data$Sample.Num, y = .data$Z.score, label = .data$Sample.Name)
  ) +
    ggplot2::geom_text(size = 4, colour = "red") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = yinter)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::xlab("Sample Number") +
    ggplot2::ylab("Z score") +
    ggplot2::ggtitle("Sample Connectivity") +
    design_mytheme(axis_angle = 0, axis_text_size = 12, axis_title_size = 2)

  if (show_plot) print(p)

  if (save) {
    ggplot2::ggsave(
      filename = paste0(index, "-2-connectivityplot.pdf"),
      plot = p,
      width = 8, height = 8,
      path = path$folder_name
    )
  }

  names_eset_rmout <- colnames(eset)[abs(connectivity.zscore) > abs(yinter)]

  cli::cli_alert_info("When yinter = {yinter}")
  cli::cli_alert_info("Potential outliers: {.val {names_eset_rmout}}")

  names_eset_rmout
}
