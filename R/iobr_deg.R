





#' Differential expressed analyses
#'
#' @param eset Expression set object containing the matrix of gene expression data. Default is NULL.
#' @param annoation The annotation to be used. Default is annotation_grch38.
#' @param pdata The data frame containing the sample information.
#' @param group_id Specifies the ID of the group. Default is "group".
#' @param pdata_id Specifies the ID of the pdata. Default is "ID".
#' @param array Logical value to specify whether the 'eset' is a microarray data or not. Default is TRUE.
#' @param method The method used for DEG analysis: either "DESeq2" or "limma". Default is "DESeq2".
#' @param contrast Specifies groups to apply DESeq2 results() function, it is typically applied to groups of interest in terms of an outcome variable of interest. Default is c("deg_group","High","Low").
#' @param padj_cutoff Specifies the cutoff of adjusted p-values. Default is 0.01.
#' @param logfc_cutoff Specifies the cutoff of log2 Fold Change. Default is 0.5.
#' @param volcano_plot Logical value to specify whether to plot volcano plots or not. Default is FALSE.
#' @param col_volcano Specifies the color of the volcano plots. Default is 1.
#' @param heatmap Logical value to specify whether to plot heatmaps or not. Default is TRUE.
#' @param col_heatmap Specifies the color of the heatmaps. Default is 1.
#' @param path
#'
#' @importFrom limma lmFit
#' @author Dongqiang Zeng
#' @return DEG object
#' @export
#'
#' @examples
#'
iobr_deg<-function(eset         ,
                   annoation    = annotation_grch38,
                   pdata,
                   group_id    = "group",
                   pdata_id     = "ID",
                   array        = TRUE,
                   method       = "DESeq2",
                   contrast     = c("deg_group","High","Low"),
                   path         = NULL,
                   padj_cutoff  = 0.01,
                   logfc_cutoff = 0.5,
                   volcano_plot = FALSE,
                   col_volcano  = 1,
                   heatmap      = TRUE,
                   col_heatmap  = 1){




  ############################################
  if(is.null(path)){
    # set path to store enrichment analyses result
    path <-creat_folder(paste0("1-result-of-DEGs"))
  }
  abspath<-path$abspath
  ############################################

  colnames(pdata)[which(colnames(pdata)==pdata_id)]<-"ID"
  colnames(pdata)[which(colnames(pdata)==group_id)]<-"deg_group"

  if(group_id=="group3"){
    pdata<-pdata[!pdata$deg_group=="Middle",]
  }
  pdata<-as.data.frame(pdata)
  pdata<-pdata[!is.na(pdata$deg_group),]
  pdata<-pdata[!pdata$deg_group=="NA",]
  pdata<-pdata[pdata$ID%in%colnames(eset),]
  eset<-eset[, colnames(eset)%in%pdata$ID]
  pdata<-pdata[match(colnames(eset),pdata$ID),]
  ########################################

  if(array) eset <- preprocessCore::normalize.quantiles(as.matrix(eset), keep.names = TRUE)

  if(method == "DESeq2"){
    eset<-round(eset,0)
    dds <-DESeq2::DESeqDataSetFromMatrix(countData = eset,
                                         colData = pdata,
                                         design= ~ deg_group)

    #过滤一些low count的数据,起码五分之一的样本有表达
    dds <- dds[rowSums(counts(dds)) > ncol(eset)/5, ]
    dds <- DESeq(dds) #,parallel = T
    # contrast<-c("deg_group","High","Low")
    res_tidy <- results(dds,tidy = T,contrast = contrast)
    res_tidy <- as.data.frame(res_tidy)
    # print(summary(res_tidy))
    #####################################
    message(paste0("Counts of gene: Adj.pvalue < 0.001: >>>  ",sum(res_tidy$padj < 0.001, na.rm=TRUE)))
    message(paste0("Counts of gene: Adj.pvalue < 0.05: >>>  ",sum(res_tidy$padj < 0.05, na.rm=TRUE)))
    message(paste0("Counts of gene: Adj.pvalue < 0.1: >>>  ",sum(res_tidy$padj < 0.1, na.rm=TRUE)))
    message(paste0("Counts of gene: Adj.pvalue < 0.25: >>>  ",sum(res_tidy$padj < 0.25, na.rm=TRUE)))


    DEG<-res_tidy[order(res_tidy$padj,decreasing = F),]
    print(head(DEG))

    if(!is.null(annotation)){
      # DEG$row<-substring(DEG$row,1,15)
      DEG<-merge(DEG, annoation, by.x="row",by.y="ensgene",all = FALSE)
    }


    DEG$sigORnot <- ifelse(DEG$log2FoldChange > logfc_cutoff & DEG$padj < padj_cutoff, "Up_regulated",
                           ifelse(DEG$log2FoldChange < -logfc_cutoff & DEG$padj < padj_cutoff, "Down_regulated", "NOT"))
    DEG$label <- ifelse(abs(DEG$log2FoldChange) > logfc_cutoff & DEG$padj < padj_cutoff, "Both",
                        ifelse(DEG$padj < padj_cutoff, "Significant",
                               ifelse(abs(DEG$log2FoldChange) >= logfc_cutoff, paste0("log2FC >= ",logfc_cutoff), "NOT")))
    #######################################################
    print(head(DEG))

    #######################################################
    contrast<-contrast[2:3]
    # contrast<-contrast[order(contrast)]
    level1<-as.character(contrast[1])
    level2<-as.character(contrast[2])
    aa<-as.character(pdata[pdata$deg_group == level1,"ID"])
    aa<-aa[aa%in%colnames(eset)]
    #######################################################
    bb<-as.character(pdata[pdata$deg_group == level2,"ID"])
    bb<-bb[bb%in%colnames(eset)]
    #######################################################
    message(paste(">>>>"))
    message(paste0("ID of ",level1," group: "))
    message(paste(aa,collapse = ", "))
    message(paste(">>>>"))
    message(paste0("ID of ",level2," group: "))
    message(paste(bb,collapse = ", "))
    #########################################################
    eset2 <- normalizeBetweenArrays(eset,method = "quantile")
    eset2<-tibble::rownames_to_column(as.data.frame(eset2),var = "ID")
    # eset2$ID<-substring(eset2$ID,1,15)
    ####################################################
    meancounts<-eset2 %>%
      mutate(mean_h=(rowSums(.[,aa]))/length(aa))%>%
      mutate(mean_l=(rowSums(.[,bb]))/length(bb))%>%
      dplyr:: select(ID,mean_h,mean_l)


    # if(!is.null(annotation)) meancounts$ID <-substring(meancounts$ID,1,15)
    ##################################################


    DEG<-merge(DEG,meancounts,by.x="row",by.y="ID",all = FALSE)
    DEG<-tibble::as_tibble(DEG)
    DEG<-DEG[order(DEG$padj,decreasing = F),]
    print(head(DEG))

  }

  if(method == "limma"){
    library(limma)
    pdata$deg_group<-ifelse(pdata$deg_group==contrast[2],"group1","group2")
    ################################
    message(paste0("group1 = ", contrast[2]))
    message(paste0("group2 = ",  contrast[3]))
    ################################
    design <- model.matrix(~0+factor(pdata$deg_group))
    colnames(design)=levels(factor(pdata$deg_group))
    rownames(design)=colnames(eset)
    ################################
    contrast.matrix<- makeContrasts(group1 - group2, levels = design)
    fit <- lmFit(eset, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    DEG <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

    DEG<-rownames_to_column(DEG,var = "symbol")
    colnames(DEG)[which(colnames(DEG)=="logFC")]<-"log2FoldChange"
    colnames(DEG)[which(colnames(DEG)=="adj.P.Val")]<-"padj"
    colnames(DEG)[which(colnames(DEG)=="P.Value")]<-"pvalue"
    #################################
    DEG<- tibble::as_tibble(DEG)
    DEG$sigORnot <- ifelse(DEG$log2FoldChange > logfc_cutoff & DEG$padj < padj_cutoff, "Up_regulated",
                           ifelse(DEG$log2FoldChange < -logfc_cutoff & DEG$padj < padj_cutoff, "Down_regulated", "NOT"))
    DEG$label <- ifelse(abs(DEG$log2FoldChange) > logfc_cutoff & DEG$padj < padj_cutoff, "Both",
                        ifelse(DEG$padj < padj_cutoff, "Significant",
                               ifelse(abs(DEG$log2FoldChange) >= logfc_cutoff, paste0("log2FC >= ",logfc_cutoff), "NOT")))

    #######################################################

    # contrast<-contrast[2:3]
    # # contrast<-contrast[order(contrast)]
    # level1<-as.character(contrast[1])
    # level2<-as.character(contrast[2])
    aa<-as.character(pdata[pdata$deg_group == "group1","ID"])
    aa<-aa[aa%in%colnames(eset)]
    #######################################################
    bb<-as.character(pdata[pdata$deg_group == "group2","ID"])
    bb<-bb[bb%in%colnames(eset)]
    ########################################################
    #########################################################
    eset2<-rownames_to_column(as.data.frame(eset),var = "ID")
    ####################################################
    meancounts<-eset2%>%mutate(mean_group1=(rowSums(.[,aa]))/length(aa))%>%
      mutate(mean_group2=(rowSums(.[,bb]))/length(bb))%>%
      dplyr:: select(ID,mean_group1,mean_group2)
    DEG<-merge(DEG, meancounts, by.x="symbol",by.y="ID",all = FALSE)
    # DEG<-tbl_df(DEG)
    DEG<-DEG[order(DEG$padj,decreasing = F),]
    DEG<-tibble::as_tibble(DEG)

    colnames(DEG)[which(colnames(DEG)=="mean_group1")] <- contrast[2]
    colnames(DEG)[which(colnames(DEG)=="mean_group2")] <- contrast[3]
    print(head(DEG))

  }

  save(DEG,file = paste0(abspath, "1-DEGs.RData"))
  writexl::write_xlsx(DEG, paste0(abspath,"2-DEGs.xlsx"))
  return(DEG)

}
