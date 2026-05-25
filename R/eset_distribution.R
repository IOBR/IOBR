#' Visualize Expression Set Distribution
#'
#' @description
#' Generates boxplots and density plots to analyze the distribution of
#' expression values in an expression set. Useful for quality control and
#' assessing data normalization.
#'
#' @param eset Expression matrix or data frame with genes in rows and samples
#'   in columns.
#' @param quantile Integer specifying the divisor for sampling columns.
#'   Default is 3 (samples 1/3 of columns).
#' @param log Logical indicating whether to perform log2 transformation.
#'   Default is `TRUE`.
#' @param project Optional output directory path for saving files. If `NULL`,
#'   no files are saved. Default is `NULL`.
#'
#' @return Invisibly returns `NULL`. If `project` is provided, saves PNG files to disk.
#'
#' @export
#'
#' @examples
#' \donttest{
#' eset_stad <- load_data("eset_stad")
#' anno_rnaseq <- load_data("anno_rnaseq")
#' eset <- anno_eset(eset = eset_stad, annotation = anno_rnaseq)
#' eset_distribution(eset)
#' eset_distribution(eset, project = file.path(tempdir(), "ESET"))
#' }
eset_distribution <- function(eset, quantile = 3, log = TRUE, project = NULL) {
  if (!is.matrix(eset) && !is.data.frame(eset)) {
    cli::cli_abort("{.arg eset} must be a matrix or data frame")
  }
  if (nrow(eset) == 0 || ncol(eset) == 0) {
    cli::cli_abort("{.arg eset} must have at least one row and one column")
  }
  if (!is.numeric(quantile) || quantile <= 0) {
    cli::cli_abort("{.arg quantile} must be a positive number")
  }

  feas <- feature_manipulation(data = eset, feature = rownames(eset), is_matrix = TRUE)
  eset <- eset[rownames(eset) %in% feas, , drop = FALSE]

  n_samples <- ncol(eset)
  n_select <- max(1, round(n_samples / quantile))
  index <- sample(seq_len(n_samples), min(n_select, n_samples))
  eset1 <- eset[, index, drop = FALSE]

  if (log) {
    eset1 <- log2eset(eset1)
  }

  eset1 <- t(eset1)
  eset1 <- as.data.frame(eset1)
  eset1$Sample.Name <- rownames(eset1)

  rlang::check_installed("reshape2")
  eset.melt <- reshape2::melt(eset1, id.vars = "Sample.Name")
  colnames(eset.melt)[2:3] <- c("Symbol", "Intensity")

  if (!is.null(project)) {
    path <- creat_folder(project)
  }

  p <- ggplot2::ggplot(eset.melt, ggplot2::aes(x = .data$Sample.Name, y = .data$Intensity)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 23),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1.0)
    ) +
    ggplot2::ggtitle(
      "Normalized signal intensity",
      paste0(
        "Patient No. is: ", n_samples, "; ",
        "Number of features is: ", nrow(eset), ";  ",
        "Maximum is: ", round(max(eset, na.rm = TRUE), 2)
      )
    ) +
    ggplot2::ylab("Intensity") +
    ggplot2::xlab("Sample") +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_text(size = 18, hjust = 0.1, face = "italic", color = "black")
    )

  p2 <- ggplot2::ggplot(eset.melt, ggplot2::aes(.data$Intensity, group = .data$Sample.Name)) +
    ggplot2::geom_density() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle("Histogram of Log2 Expression") +
    ggplot2::ylab("Density") +
    ggplot2::xlab("Log2 Expression")

  if (!is.null(project)) {
    ggplot2::ggsave(
      filename = paste0("1-", project, "-boxplot.png"),
      plot = p,
      width = 15, height = 8,
      path = path$folder_name
    )

    ggplot2::ggsave(
      filename = paste0("2-", project, "-Densityplot.png"),
      plot = p2,
      width = 9, height = 6,
      path = path$folder_name
    )
  }

  invisible(NULL)
}
