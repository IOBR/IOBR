#' Enrichment Bar Plot with Two Directions
#'
#' Creates a bar plot visualizing enrichment results for up-regulated and down-regulated terms,
#' using -log10(p-values) to indicate significance.
#'
#' @param up_terms Data frame for up-regulated terms.
#' @param down_terms Data frame for down-regulated terms.
#' @param title Plot title. Default is "Gene Ontology Enrichment".
#' @param width_wrap Maximum width for wrapping pathway names. Default is 30.
#' @param palette Color palette. Default is "jama".
#' @param terms Column name for term descriptions. Default is "Description".
#' @param pvalue Column name for p-values. Default is "pvalue".
#' @param group Column name for group indicator. Default is "group".
#' @param font_terms Font size for axis labels. Default is 15.
#'
#' @return A ggplot object of the enrichment bar plot.
#' @export
#' @author Dongqiang Zeng
#' @examples
#' up_terms <- data.frame(Description = c("Pathway1", "Pathway2"), pvalue = c(0.001, 0.01))
#' down_terms <- data.frame(Description = c("Pathway4", "Pathway5"), pvalue = c(0.005, 0.02))
#' p <- enrichment_barplot(up_terms = up_terms, down_terms = down_terms,
#'                         title = "Custom Enrichment Plot")
enrichment_barplot <- function(up_terms, down_terms, terms = "Description", pvalue = "pvalue", group = "group", palette = "jama", title = "Gene Ontology Enrichment", width_wrap = 30, font_terms = 15) {
  up_terms <- as.data.frame(up_terms)
  down_terms <- as.data.frame(down_terms)
  dat <- rbind(as.data.frame(up_terms), as.data.frame(down_terms))
  ############################################
  colnames(dat)[which(colnames(dat) == terms)] <- "terms"
  colnames(dat)[which(colnames(dat) == pvalue)] <- "pvalue"

  if (!group %in% colnames(dat)) {
    dat$group <- ifelse(dat$terms %in% up_terms[, terms], 1, -1)
  } else {
    colnames(dat)[which(colnames(dat) == group)] <- "group"
  }

  dat$terms <- gsub(dat$terms, pattern = "\\_", replacement = " ")
  dat$pvalue <- -log10(as.numeric(dat$pvalue))
  dat$pvalue <- dat$pvalue * dat$group
  dat$pvalue <- as.numeric(dat$pvalue)
  dat <- dat[order(dat$pvalue, decreasing = F), ]

  # print(dat)
  color <- palettes(category = "box", palette = palette, alpha = 1, show_col = FALSE)

  p <- ggplot(dat, aes(x = reorder(terms, order(pvalue, decreasing = F)), y = pvalue, fill = group)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = color[1], high = color[2], guide = FALSE) +
    scale_x_discrete(name = "Pathway names", labels = function(x) str_wrap(x, width = width_wrap)) +
    scale_y_continuous(name = "log10(P.value)") +
    coord_flip() +
    theme_light() +
    theme(
      plot.title = element_text(hjust = 0),
      axis.title = element_text(size = rel(1.3)),
      axis.text.x = element_text(face = "plain", size = font_terms, angle = 0, color = "black"),
      axis.text.y = element_text(face = "plain", size = font_terms, angle = 0, color = "black")
    ) +
    ggtitle(title)

  return(p)
}
