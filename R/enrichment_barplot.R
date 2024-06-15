





#' Enrichment bar plot with two directions
#'
#' This function creates a bar plot to visualize enrichment results, displaying both up-regulated and down-regulated terms or genes. The bar plot highlights the significance of each term or gene using -log10(p-value).
#'
#' @param up_terms A data frame containing the information for up-regulated terms
#' @param down_terms A data frame containing the information for down-regulated terms or genes.
#' @param title The title for the barplot. By default, it is set to "Gene Ontology Enrichment".
#' @param width_wrap The maximum width for wrapping the pathway names on the x-axis labels. By default, it is set to 30 characters.
#' @param palette The color palette to use for the barplot. By default, it is set to "jama".
#' @param terms The column name in the data frames that specifies the term or gene description. By default, it is set to "Description".
#' @param pvalue The column name in the data frames that specifies the p-value. By default, it is set to "pvalue".
#' @param group The column name in the data frames that specifies the group or condition. By default, it is set to "group". If this column is not present in the data frames, the function will automatically assign a group based on the up_terms and down_terms.  If the group is "Upregulated", it is assigned a value of 1, otherwise it is assigned a value of -1.
#' @param font_terms The font size for the x-axis and y-axis labels. By default, it is set to 15.
#'
#' @author Dongqiang Zeng
#' @return A ggplot object representing the enrichment bar plot.
#' @export
#'
#' @examples
#' # Example data frames
#' up_terms <- data.frame(
#'   Description = c("Pathway1", "Pathway2", "Pathway3"),
#'   pvalue      = c(0.001, 0.01, 0.05))
#' 
#' down_terms <- data.frame(
#'   Description = c("Pathway4", "Pathway5", "Pathway6"),
#'   pvalue      = c(0.005, 0.02, 0.03))
#' 
#' # Create enrichment bar plot using the function
#' p <- enrichment_barplot(
#'   up_terms   = up_terms,
#'   down_terms = down_terms,
#'   title      = "Custom Enrichment Plot")
#' 
#' print(p)

enrichment_barplot <- function(up_terms, down_terms, terms = "Description", pvalue = "pvalue", group = "group", palette = "jama", title = "Gene Ontology Enrichment", width_wrap = 30, font_terms = 15){

  up_terms<- as.data.frame(up_terms)
  down_terms<-as.data.frame(down_terms)
  dat<-rbind(as.data.frame(up_terms), as.data.frame(down_terms))
  ############################################
  colnames(dat)[which(colnames(dat)== terms )]<- "terms"
  colnames(dat)[which(colnames(dat)== pvalue )]<- "pvalue"

  if(!group%in%colnames(dat)){
    dat$group <-ifelse( dat$terms %in% up_terms[,terms], 1, -1)
  }else{
    colnames(dat)[which(colnames(dat)== group )]<- "group"
  }

  dat$terms<- gsub(dat$terms, pattern = "\\_",replacement = " ")
  dat$pvalue <- -log10(as.numeric(dat$pvalue))
  dat$pvalue<- dat$pvalue*dat$group
  dat$pvalue<- as.numeric(dat$pvalue)
  dat<-dat[order(dat$pvalue,decreasing = F),]

  # print(dat)
  color<- palettes(category = "box", palette = palette, alpha = 1, show_col = FALSE)

  p<- ggplot(dat, aes(x=reorder(terms,order(pvalue, decreasing = F)), y=pvalue, fill=group)) +
    geom_bar(stat="identity") +
    scale_fill_gradient(low = color[1],high = color[2],guide = FALSE) +
    scale_x_discrete(name ="Pathway names",labels=function(x) str_wrap(x, width = width_wrap)) +
    scale_y_continuous(name ="log10(P.value)") +
    coord_flip() + theme_light()+
    theme(plot.title = element_text(hjust = 0),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x= element_text(face="plain",size=font_terms, angle=0,color="black"),
          axis.text.y= element_text(face="plain",size=font_terms, angle=0,color="black"))+
    ggtitle(title)

  return(p)

}


