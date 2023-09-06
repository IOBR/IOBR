


#' sig_gsea - Perform Gene Set Enrichment Analysis
#'
#' @description The sig_gsea function conducts Gene Set Enrichment Analysis (GSEA) to identify significant gene sets based on differential gene expression data. It allows for customization of various parameters to tailor the analysis to specific needs. This function performs GSEA using the fgsea package and provides visualizations and results in the form of tables and plots. It also supports the utilization of user-defined gene sets or the use of predefined gene sets from the Molecular Signatures Database (MSigDB). The function further allows for customization of parameters such as organism, gene symbol type, visualization color palette, and significance thresholds. The results can be saved in Excel format, and plots can be saved in various image formats.
#' @param deg Differential expressed genes object, typically a data frame that contains gene symbols, log fold changes, and other relevant information.
#' @param genesets This parameter allows you to specify a custom set of gene sets to be used in the enrichment analysis. If not provided, the function will use the gene sets available in the "msigdb" database based on the selected organism.
#' @param path The path parameter represents the location where the enrichment analysis results will be stored. If not specified, a default path named "1-GSEA-result" will be created in the current working directory.
#' @param gene_symbol his parameter specifies the column name in the deg data frame that contains the gene symbols. The default value is "symbol".
#' @param logfc Specifies the column name in the deg data frame that contains the log fold change values. The default value is "log2FoldChange".
#' @param org This parameter is used to select the organism for which the enrichment analysis will be performed. The options are "hsa" for Homo sapiens and "mus" for Mus musculus.
#' @param msigdb A logical parameter indicating whether to use the gene sets from the "msigdb" database. If set to TRUE, the function will retrieve gene sets from "msigdb" based on the selected organism and category.
#' @param category Specifies the category of gene sets to be used from the "msigdb" database. The default category is "H", representing Hallmark gene sets.
#' @param subcategory Allows you to specify a subcategory of gene sets from the "msigdb" database. If not provided, all gene sets within the selected category will be used.
#' @param palette_bar Specifies the color palette for the barplot used to visualize the enriched gene sets. The default value is "nrc".
#' @param palette_gsea Specifies the color palette for the GSEA plots. The default value is 2.
#' @param show_bar Specifies the number of enriched gene sets to show in the barplot. The default value is 10.
#' @param show_col A logical parameter indicating whether to show the color names in the barplot. The default value is FALSE.
#' @param show_plot A logical parameter indicating whether to display the GSEA plots. The default value is FALSE.
#' @param show_gsea Specifies the number of most significant gene sets to show in the GSEA plots. The default value is 8.
#' @param plot_single_sig A logical parameter indicating whether to plot each significant gene set separately. The default value is TRUE.
#' @param project Specifies the name of the project or category for the analysis. If not provided, it will be set as "custom_sig".
#' @param minGSSize Specifies the minimum gene set size to consider for enrichment analysis. Gene sets below this size will be excluded. The default value is 10.
#' @param maxGSSize Specifies the maximum gene set size to consider for enrichment analysis. Gene sets above this size will be excluded. The default value is 500.
#' @param verbose A logical parameter indicating whether to display additional information and messages during the analysis. The default value is TRUE.
#' @param seed A logical parameter indicating whether to use a random seed for reproducibility in the analysis. The default value is FALSE.
#' @param fig.type Specifies the file type for saving the GSEA plots. The default value is "pdf".
#' @param show_path_n Specifies the number of paths to show in the GSEA plots. The default value is 20.
#' @param print_bar Default is TRUE
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
sig_gsea <-function(deg,
                    genesets          = NULL,
                    path              = NULL,
                    gene_symbol       = "symbol",
                    logfc             = "log2FoldChange",
                    org               = "hsa",
                    msigdb            = TRUE,
                    category          = "H",
                    subcategory       = NULL,
                    palette_bar       = "jama",
                    palette_gsea      = 2,
                    show_bar          = 10,
                    show_col          = FALSE,
                    show_plot         = FALSE,
                    show_gsea         = 8,
                    show_path_n       = 20,
                    plot_single_sig   = FALSE,
                    project           = "custom_sig",
                    minGSSize         = 10,
                    maxGSSize         = 500,
                    verbose           = TRUE,
                    seed              = FALSE,
                    fig.type          = "pdf",
                    print_bar         = TRUE){


  # set path to store enrichment analyses result
  if(is.null(path)){
    file_store<-paste0("1-GSEA-result")
  }else{
    file_store<-path
  }

  if ( ! file.exists(file_store) ) dir.create(file_store)
  abspath<-paste(getwd(),"/",file_store,"/",sep ="" )
  #################################################

  # if(!is.null(input)) (load(paste0(file_source,"/",input)))

  if(gene_symbol!="symbol"){
    colnames(deg)[which(colnames(deg)=="symbol")]<-"symbol_subs"
  }

  colnames(deg)[which(colnames(deg)==gene_symbol)]<-"symbol"
  colnames(deg)[which(colnames(deg)==logfc)]<-"logfc"
  ##################################

  message("`>>>--- Parametar org must be one of: hsa or mus`")
  if(org=="hsa"){

    database<-"org.Hs.eg.db"
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
    entrizid <-clusterProfiler:: bitr(deg$symbol, fromType = "SYMBOL",
                     toType = c("GENENAME", "ENTREZID"),
                     OrgDb = database)
  }else if(org =="mus"){
    database<-"org.Mm.eg.db"
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
    entrizid <-clusterProfiler:: bitr(deg$symbol, fromType = "SYMBOL",
                     toType = c("GENENAME", "ENTREZID"),
                     OrgDb = database)
  }

  #################################
  deg<-merge(deg,entrizid,by.x="symbol",by.y="SYMBOL",all.x=T,all.y=F)
  #还是发现通过org.hs.eg.db注释所获得的ID更多
  ##################################
  gene_id_logfc <- deg %>%  dplyr::select(ENTREZID, logfc) %>%
    dplyr::distinct(ENTREZID,.keep_all = T) %>%
    filter(!is.na(ENTREZID)) %>%
    filter(!is.na(logfc)) %>%
    arrange(desc(logfc))
  #按FC降序
  genelist<-gene_id_logfc$logfc
  names(genelist)<-gene_id_logfc$ENTREZID
  genelist<-genelist[order(genelist,decreasing = T)]
  #########################################


  if(is.null(genesets)){

    if(org=="hsa"){
      species<- "Homo sapiens"
    }else if(org=="mus"){
      species<- "Mus musculus"
    }
    ##################################################
    message(">>>---- Categories that can be choosed... ")

    m_df = msigdbr::msigdbr(species = species)
    a <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
    print(as.data.frame(a))

    ##################################################
    if(is.null(category)){
      message(">>>--- Category is NULL, default is Hallmark gene sets...")
      category = "H"
    }
    term2genes <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
    term2genes <- term2genes %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

  }else{

    # if(org=="mus"){
    #   genesets<-lapply(genesets, function(x) mus2human_gene)
    # }

    term2genes<- output_sig(genesets,file.name = "sig")
    file.remove("sig.csv")
    term2genes<-reshape2:: melt(as.matrix(term2genes))
    term2genes<-term2genes[,-1]
    colnames(term2genes)<-c("gs_name","symbol")
    category<-project


    if(org=="hsa"){

      database<-"org.Hs.eg.db"
      entrizid <- bitr(term2genes$symbol, fromType = "SYMBOL",
                       toType = c("GENENAME", "ENTREZID"),
                       OrgDb = database)
    }else if(org =="mus"){
      database<-"org.Mm.eg.db"
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
      entrizid <- bitr(term2genes$symbol, fromType = "SYMBOL",
                       toType = c("GENENAME", "ENTREZID"),
                       OrgDb = database)
    }

    #################################
    term2genes<-merge(term2genes,entrizid,by.x="symbol",by.y="SYMBOL",all.x=T,all.y=F)

    term2genes<-term2genes[,c("gs_name","ENTREZID")]
    colnames(term2genes)<-c("gs_name","entrez_gene")
  }

  ####################################################
  hall_gsea <-clusterProfiler::GSEA(genelist,
                                    exponent      = 1,
                                    # nPerm         = 1000,
                                    eps           = 0,
                                    minGSSize     = minGSSize,
                                    maxGSSize     = maxGSSize,
                                    pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH",
                                    TERM2GENE     = term2genes,
                                    TERM2NAME     = NA,
                                    verbose       = verbose,
                                    seed          = FALSE,
                                    by            = "fgsea")

  if(org=="hsa"){
    hall_gsea<- DOSE:: setReadable(hall_gsea, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  }else if(org=="mus"){
    hall_gsea<- DOSE:: setReadable(hall_gsea, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
  }

  writexl::write_xlsx(as.data.frame(hall_gsea),paste0(abspath,"1-",category,"_GSEA_significant_results.xlsx"))
  # save(hall_gsea,file = paste(abspath,"2-",category,"_GSEA_result.RData",sep = ""))

  ###########################################
  if(dim(hall_gsea)[1] > 0){

    if(dim(hall_gsea)[1]<show_gsea){
      paths<-rownames(hall_gsea@result)
    }else{
      paths<-rownames(hall_gsea[1:show_gsea,])
    }

    print(paste0(">>>--- Most significant gene sets: "))
    print(paste(1:length(paths),paths,collapse = "; "))
    ###############################################

    # gseacol<- palettes(category = "random", palette = palette_gsea, show_col = FALSE, show_message = F)
    gseacol <- get_cols(palette = palette_gsea, show_col = FALSE)
    ###############################################
    GSEAPLOT<-enrichplot:: gseaplot2(x = hall_gsea,
                                    geneSetID =  paths,
                                     subplots     =c(1,2),
                                     color        = gseacol[1:length(paths)],
                                     pvalue_table = TRUE)

    ggplot2::ggsave(filename = paste0("2-",category,"_Top_",show_gsea,"_GSEA_plot.",fig.type),
                    plot = GSEAPLOT, path = file_store,
                    width = 11, height = 7, dpi = 300)

    if(show_plot) print(GSEAPLOT)
    ###############################################

    if(plot_single_sig) {

      paths<- rownames(hall_gsea@result)

      paths2<-paths[1:show_path_n]
      paths2  <- paths2[!is.na(paths2)]

      for (i in 1:length(paths2)) {
        single_path<-paths2[i]
        message(paste0(">>>--- Processing signature: ", single_path))


        ###############################################
        gseaplot<-enrichplot:: gseaplot2(x = hall_gsea,
                                         geneSetID = single_path,
                                         subplots     =c(1,2),
                                         color        = gseacol[i],
                                         pvalue_table = TRUE)
        gseaplot<-gseaplot
        # +design_mytheme(axis_angle = 0, hjust = 0.5)

        if(show_plot) print(gseaplot)
        ggplot2::ggsave(filename = paste0(i+4,"-GSEA_plot-",single_path,".",fig.type), plot = gseaplot,
                        path = file_store, width = 11,height = 7.5, dpi = 300)

      }

    }

    #' 画出P_value的barplot
    ################################################
    down_gogo<-hall_gsea[hall_gsea$pvalue<0.05 & hall_gsea$enrichmentScore < 0,]
    if(!dim(down_gogo)[1]==0) down_gogo$group=-1
    down_gogo<-down_gogo[1:show_path_n,]
    down_gogo<-down_gogo[!is.na(down_gogo$p.adjust),]

    ################################################
    up_gogo<-hall_gsea[hall_gsea$pvalue<0.05 & hall_gsea$enrichmentScore > 0,]
    if(!dim(up_gogo)[1]==0) up_gogo$group= 1
    up_gogo<-up_gogo[1:show_path_n,]
    up_gogo<-up_gogo[!is.na(up_gogo$p.adjust),]


   if(print_bar){
    ################################################

    gsea_bar<-enrichment_barplot(up_terms   = up_gogo,
                                 down_terms = down_gogo,
                                 palette    = palette_bar,
                                 title      = "GSEA-Enrichment",
                                 width_wrap = 30)

    n_bar <- c(dim(up_gogo)[1]+ dim(down_gogo)[1])
    height_bar <- 0.5*n_bar + 3
   ggplot2:: ggsave(filename = paste0('3-', category,'_GSEA_barplot.',fig.type), plot = gsea_bar,
                    path = file_store, width = 6,height = height_bar)

    if(show_plot) print(gsea_bar)
    ######################################################
  }

    hall_gsea <- as.tibble(hall_gsea)
    #' restore data
    hall_gsea_result<-list(up = hall_gsea[hall_gsea$pvalue<0.05 & hall_gsea$enrichmentScore > 0,],
                           down = hall_gsea[hall_gsea$pvalue<0.05 & hall_gsea$enrichmentScore < 0,],
                           all = hall_gsea,
                           plot_top = GSEAPLOT)
    save(hall_gsea_result, file = paste(abspath,"0-",category,"_combined_GSEA_result.RData",sep = ""))
    #######################################################
  }else{
    hall_gsea_result<-hall_gsea
    print(paste0(">>>--- No terms enriched under specific pvalue Cutoff = 0.05"))
  }


  return(hall_gsea_result)

}
