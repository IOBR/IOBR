





#' Enrichment barplot with two directions
#'
#' @param up_terms Enrichment result 'data frame' of up regulated genes
#' @param down_terms Enrichment result 'data frame' of down regulated genes
#' @param title plot title of barplot
#' @param width_wrap wrap the name of terms if it is too long
#' @param palette
#' @param terms
#' @param pvalue
#' @param group
#' @param font_terms
#'
#' @author Dongqiang Zeng
#' @return
#' @export
#'
#' @examples
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

  print(dat)
  color<- palettes(category = "box", palette = palette, alpha = 1)

  p<- ggplot(dat, aes(x=reorder(terms,order(pvalue, decreasing = F)), y=pvalue, fill=group)) +
    geom_bar(stat="identity") +
    scale_fill_gradient(low = color[1],high = color[2],guide = FALSE) +
    scale_x_discrete(name ="Pathway names",labels=function(x) str_wrap(x, width = width_wrap)) +
    scale_y_continuous(name ="log10(P.value)") +
    coord_flip() + theme_light()+
    theme(plot.title = element_text(hjust = 0),
          axis.title=element_text(size=rel(3)),
          axis.text.x= element_text(face="plain",size=font_terms, angle=0,color="black"),
          axis.text.y= element_text(face="plain",size=font_terms, angle=0,color="black"))+
    ggtitle(title)

  return(p)

}


