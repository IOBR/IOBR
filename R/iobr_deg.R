





#' iobr_deg
#'
#' @description This function performs differential expression analysis on gene expression data using the DESeq2 or limma method. It filters low count data, calculates fold changes and adjusted p-values, and identifies differentially expressed genes (DEGs) based on specified cutoffs. It also provides optional visualization tools such as volcano plots and heatmaps. The function saves the results as an RData file and an Excel file.
#'
#' @param eset The gene expression data matrix
#' @param annotation (Optional) The annotation object for mapping gene IDs to gene names. Default is NULL.
#' @param pdata The DataFrame containing sample information and group labels.
#' @param group_id The column name in "pdata" that specifies the group labels. Default value is "group".
#' @param pdata_id The column name in "pdata" that specifies the sample IDs. Default value is "ID".
#' @param array Logical value indicating whether to normalize the gene expression data using quantile normalization. Default value is TRUE.
#' @param method The method to be used for differential expression analysis. Options are "DESeq2" or "limma". Default value is "DESeq2".
#' @param contrast A character vector specifying the contrast for the analysis.  The first and second elements are the names of the two groups being compared. Default value is  "High" and "Low".
#' @param padj_cutoff The adjusted p-value cutoff for determining significance. Default value is 0.01.
#' @param logfc_cutoff The log2 fold change cutoff for determining significance. Default value is 0.5.
#' @param volcano_plot Logical value indicating whether to generate a volcano plot. Default value is FALSE.
#' @param col_volcano The color index for plotting the volcano plot. Default value is 1.
#' @param heatmap Logical value indicating whether to generate a heatmap. Default value is TRUE.
#' @param col_heatmap The color index for plotting the heatmap. Default value is 1.
#' @param path path to save results
#' @param parallel default is FALSE
#' @param id_anno Identifier of annotation data matched to row names of eset
#'
#' @importFrom limma lmFit
#' @author Dongqiang Zeng
#' @return DEG object
#' @export
#'
#' @examples
#'
iobr_deg<-function(eset,
                   annotation    = NULL,
                   id_anno       = NULL,
                   pdata,
                   group_id     = "group",
                   pdata_id     = "ID",
                   array        = FALSE,
                   method       = "DESeq2",
                   contrast     = c("High","Low"),
                   path         = NULL,
                   padj_cutoff  = 0.01,
                   logfc_cutoff = 0.5,
                   volcano_plot = FALSE,
                   col_volcano  = 1,
                   heatmap      = TRUE,
                   col_heatmap  = 1,
                   parallel     = FALSE){

  ############################################
  if(is.null(path)){
    # set path to store enrichment analyses result
    path <-creat_folder(paste0("Result-of-DEGs"))
  }else{
    path <- creat_folder(path)
  }
  abspath<-path$abspath
  ############################################

  colnames(pdata)[which(colnames(pdata)==pdata_id)]<-"ID"
  colnames(pdata)[which(colnames(pdata)==group_id)]<-"deg_group"

  message(">>>== Matching grouping information and expression matrix \n")
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

    message(">>>== DEGseq2 (method) was selected for differential gene analysis of RNAseq \n")

    message(">>>== Please ensure that `eset` is a count expression matrix \n")
    #########################################
    eset<-round(eset,0)
    dds <-DESeq2::DESeqDataSetFromMatrix(countData = eset,
                                         colData = pdata,
                                         design= ~ deg_group)

    #过滤一些low count的数据,起码五分之一的样本有表达
    dds <- dds[rowSums(counts(dds)) > ncol(eset)/5, ]
    dds <- DESeq(dds, parallel = parallel) #,parallel = T


    contrast<-c("deg_group", contrast)
    res_tidy <- results(dds,tidy = T,contrast = contrast)
    res_tidy <- as.data.frame(res_tidy)
    # print(summary(res_tidy))
    #####################################
    message(">>>== Summary of differential gene analysis results \n")
    message(paste0("Counts of gene: Adj.pvalue < 0.001:  ",sum(res_tidy$padj < 0.001, na.rm=TRUE)))
    message(paste0("Counts of gene: Adj.pvalue < 0.05:  ",sum(res_tidy$padj < 0.05, na.rm=TRUE)))
    message(paste0("Counts of gene: Adj.pvalue < 0.1:  ",sum(res_tidy$padj < 0.1, na.rm=TRUE)))
    message(paste0("Counts of gene: Adj.pvalue < 0.25:  ",sum(res_tidy$padj < 0.25, na.rm=TRUE)))

    DEG <- res_tidy[order(res_tidy$padj, decreasing = F), ]


    if(!is.null(annotation)){
      # DEG$row<-substring(DEG$row,1,15)
      colnames(annotation)[which(colnames(annotation)==id_anno)] <-"id"
      DEG <- merge(DEG, annotation, by.x="row", by.y="id",all = FALSE)

    }else{
      # print(head(DEG))
      # DEG <- rownames_to_column(DEG, var = "row")

      message(">>>== IOBR provides annotation files (`anno_grch38`) to help you annotate the results of `iobr_deg` \n")
      data("anno_grch38", package = "IOBR")
      DEG <- merge(DEG, anno_grch38, by.x="row", by.y="id",all = FALSE)

    }

    # print(head(DEG))
    DEG$sigORnot <- ifelse(DEG$log2FoldChange > logfc_cutoff & DEG$padj < padj_cutoff, "Up_regulated",
                           ifelse(DEG$log2FoldChange < -logfc_cutoff & DEG$padj < padj_cutoff, "Down_regulated", "NOT"))
    DEG$label <- ifelse(abs(DEG$log2FoldChange) > logfc_cutoff & DEG$padj < padj_cutoff, "Both",
                        ifelse(DEG$padj < padj_cutoff, "Significant",
                               ifelse(abs(DEG$log2FoldChange) >= logfc_cutoff, paste0("log2FC >= ",logfc_cutoff), "NOT")))
    #######################################################
    # print(head(DEG))
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
    if(length(aa)<= 20 & length(bb) <= 20) {
      message(paste(">>>>"))
      message(paste0("ID of ",level1," group: "))
      message(paste(aa, collapse = ", "))
      message(paste(">>>>"))
      message(paste0("ID of ",level2," group: "))
      message(paste(bb, collapse = ", "))
    }
    #########################################################
    eset2 <- preprocessCore::normalize.quantiles(as.matrix(eset), keep.names = TRUE)
    eset2<-tibble::rownames_to_column(as.data.frame(eset2), var = "ID")
    # eset2$ID<-substring(eset2$ID,1,15)
    ####################################################
    # meancounts<-eset2 %>%
    #   mutate(mean_h=(rowSums(.[,aa]))/length(aa))%>%
    #   mutate(mean_l=(rowSums(.[,bb]))/length(bb))%>%
    #   dplyr:: select(ID,mean_h, mean_l)
    message(paste0("group1 = ", contrast[1]))
    message(paste0("group2 = ",  contrast[2]))

    meancounts<-eset2%>%mutate(mean_group1=(rowSums(.[,aa]))/length(aa))%>%
      mutate(mean_group2=(rowSums(.[,bb]))/length(bb))%>%
      dplyr:: select(ID, mean_group1,mean_group2)

    # if(!is.null(annotation)) meancounts$ID <-substring(meancounts$ID,1,15)
    ##################################################
    # print(meancounts)

    DEG<-merge(DEG, meancounts, by.x="row", by.y="ID",all = FALSE)
    colnames(DEG)[which(colnames(DEG)=="mean_group1")] <- contrast[1]
    colnames(DEG)[which(colnames(DEG)=="mean_group2")] <- contrast[2]

    DEG<-tibble::as_tibble(DEG)
    DEG<-DEG[order(DEG$padj,decreasing = F),]
    print(head(DEG))

  }

  if(method == "limma"){

    message(">>>== limma was selected for differential gene analysis of Array data \n")

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
    DEG<-merge(DEG, meancounts, by.x= "symbol", by.y="ID", all = FALSE)
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
