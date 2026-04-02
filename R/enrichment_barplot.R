#' Enrichment Bar Plot with Two Directions
#'
#' @description
#' Creates a bar plot visualizing enrichment results for up-regulated and
#' down-regulated terms, using -log10(p-values) to indicate significance.
#'
#' @param up_terms Data frame for up-regulated terms.
#' @param down_terms Data frame for down-regulated terms.
#' @param terms Column name for term descriptions. Default is "Description".
#' @param pvalue Column name for p-values. Default is "pvalue".
#' @param group Column name for group indicator. Default is "group".
#' @param palette Color palette. Default is "jama".
#' @param cols Character vector. Custom colors for bars. If NULL, uses palette.
#'   Default is NULL.
#' @param title Plot title. Default is "Gene Ontology Enrichment".
#' @param width_wrap Maximum width for wrapping pathway names. Default is 30.
#' @param font_terms Font size for axis labels. Default is 15.
#'
#' @return A ggplot object of the enrichment bar plot.
#'
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' up_terms <- data.frame(
#'   Description = c("Pathway1", "Pathway2"),
#'   pvalue = c(0.001, 0.01)
#' )
#' down_terms <- data.frame(
#'   Description = c("Pathway4", "Pathway5"),
#'   pvalue = c(0.005, 0.02)
#' )
#' p <- enrichment_barplot(
#'   up_terms = up_terms,
#'   down_terms = down_terms,
#'   title = "Custom Enrichment Plot"
#' )
#' p
enrichment_barplot <- function(up_terms, down_terms,
                               terms = "Description",
                               pvalue = "pvalue",
                               group = "group",
                               palette = "jama",
                               cols = NULL,
                               title = "Gene Ontology Enrichment",
                               width_wrap = 30,
                               font_terms = 15) {
  up_terms <- as.data.frame(up_terms)
  down_terms <- as.data.frame(down_terms)
  dat <- rbind(up_terms, down_terms)

  # Rename columns for internal use
  colnames(dat)[which(colnames(dat) == terms)] <- "terms"
  colnames(dat)[which(colnames(dat) == pvalue)] <- "pvalue"

  # Set up group indicator
  if (!group %in% colnames(dat)) {
    dat$group <- ifelse(dat$terms %in% up_terms[, terms], 1, -1)
  } else {
    colnames(dat)[which(colnames(dat) == group)] <- "group"
  }

  dat$terms <- gsub(dat$terms, pattern = "_", replacement = " ")
  dat$pvalue <- -log10(as.numeric(dat$pvalue))
  dat$pvalue <- dat$pvalue * dat$group
  dat <- dat[order(dat$pvalue, decreasing = FALSE), ]

  # Get colors
  if (is.null(cols)) {
    color <- palettes(category = "box", palette = palette, show_col = FALSE)
  } else {
    if (length(cols) < 2) {
      cli::cli_abort("{.arg cols} must have at least 2 colors")
    }
    color <- cols[1:2]
  }

  ggplot2::ggplot(dat, ggplot2::aes(
    x = reorder(.data$terms, order(.data$pvalue, decreasing = FALSE)),
    y = .data$pvalue, fill = .data$group
  )) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_gradient(low = color[1], high = color[2], guide = "none") +
    ggplot2::scale_x_discrete(
      name = "Pathway names",
      labels = function(x) stringr::str_wrap(x, width = width_wrap)
    ) +
    ggplot2::scale_y_continuous(name = "-log10(P.value)") +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.3)),
      axis.text.x = ggplot2::element_text(
        face = "plain", size = font_terms, angle = 0, color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        face = "plain", size = font_terms, angle = 0, color = "black"
      )
    ) +
    ggplot2::ggtitle(title)
}
