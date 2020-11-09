


####################################
#' Integrative correlation between phenotype and features
#'
#' @param pdata_group pdata
#' @param id1 column of identifier of pdata
#' @param feature_data feature data with first column as `id2`
#' @param id2 column of identifier of feature data
#' @param target default is NULL, column of target variables if target is continuous
#' @param group column used to differentiate patient groups
#' @param is_target_continuous logical variable, if TRUE, new group will be generated based on the average or the third percentile
#' @param padj_cutoff cutoff of adjust p-value to filter features when number of features is larger than `feature_limit`
#' @param index index use to order the file name
#' @param signature_group for signature group `sig_group`, for gene_group `signature_collection` or `signature_tme`
#' @param category `signature` or `gene`
#' @param ProjectID names to define the file names
#' @param feature_limit maximal number features to plot
#' @param character_limit maximal number of feature characters
#' @param palette_box palette for box plot
#' @param palette_corplot palette for cor-plot
#' @param palette_heatmap palette for heatmap
#' @param show_heatmap_col_name logical variable, if TRUE, `tidyheatmap` will print the column name
#' @param show_col logical variable, if TRUE, color will be print and show in the R studio
#' @param show_plot logical variable, if TRUE, plots will be print
#'
#' @author Dongqiang Zeng
#' @return
#' @export
#'
#' @examples
#'
iobr_cor_plot<-function(pdata_group,
                         id1 = "ID",
                         feature_data,
                         id2 = "ID",
                         target = NULL,
                         group = "group3",
                         is_target_continuous = TRUE,
                         padj_cutoff = 1,
                         index = 1,
                         category = "signature",
                         signature_group = sig_group,
                         ProjectID = "TCGA-STAD",
                         palette_box = "nrc",
                         palette_corplot = "pheatmap",
                         palette_heatmap = 2,
                         feature_limit = 26,
                         character_limit = 30,
                         show_heatmap_col_name = FALSE,
                         show_col = FALSE,
                         show_plot = FALSE){

  if(!is.null(target)){
    file_store<-paste0(index,"-1-",ProjectID,"-",target,"-relevant-",category)
  }else{
    file_store<-paste0(index,"-1-",ProjectID,"-",group,"-relevant-",category)
  }

  if ( ! file.exists(file_store) ) dir.create(file_store)
  abspath<-paste(getwd(),"/",file_store,"/",sep ="" )

  if(is.null(names(signature_group))) signature_group<-list('signature' = signature_group)
  ##################################################
  pdata_group<-as.data.frame(pdata_group)
  colnames(pdata_group)[which(colnames(pdata_group)==id1)]<-"ID"

  if(!is.null(target)&is_target_continuous) pdata_group[,target]<-as.numeric(pdata_group[,target])

  if(is_target_continuous&!"group2"%in%colnames(pdata_group)){
    pdata_group$group2<-ifelse(pdata_group[,target]>=mean(pdata_group[,target]),"High","Low")
  }

  if(is_target_continuous&!"group3"%in%colnames(pdata_group)){
    q1<-quantile(pdata_group[,target],probs = 1/3)
    q2<-quantile(pdata_group[,target],probs = 2/3)
    pdata_group$group3<-ifelse(pdata_group[,target]<=q1,"Low",ifelse(pdata_group[,target]>=q2,"High","Middle"))
  }

  if(!is.null(target)&is_target_continuous){
    pdata_group<-pdata_group[,c("ID",target,"group2","group3")]
  }else{
    pdata_group<-pdata_group[,c("ID",group)]
  }

  # pdata_group<-pdata_group[!is.na(pdata_group[,group]),]
  # pdata_group<-pdata_group[!pdata_group[,group]=="NA",]
  #####################################################

  if(category=="gene"){
    #Gene expression data should be transposed
    feature_data<-tibble::rownames_to_column(as.data.frame(t(feature_data)),var = "ID")
  }

  feature_data<-as.data.frame(feature_data)
  colnames(feature_data)[which(colnames(feature_data)==id2)]<-"ID"
  feature_selected<-feature_manipulation(data = feature_data,feature = setdiff(colnames(feature_data),"ID"))
  feature_data<-feature_data[,colnames(feature_data)%in%c("ID",feature_selected)]

  if(!is.null(target)){
    if(target%in%colnames(feature_data)){
      feature_data<-feature_data[,-which(colnames(feature_data)==target)]
    }

  }

  if(category == "signature"){
    group_list<-signature_group
    panel<-names(signature_group)
    feature_data<-feature_data[,colnames(feature_data)%in%c("ID",unique(unlist(group_list)))]
    title.y<-"Signature score"
    title.x<-"Signatures"
  }
  if(category == "gene"){
    group_list<-signature_collection[names(signature_collection)%in%names(signature_group)]
    panel<-names(signature_group)
    feature_data<-feature_data[,colnames(feature_data)%in%c("ID",unique(unlist(group_list)))]
    title.y<-"Gene expression"
    title.x<-"Signature genes"
  }


  feature_max <- max(feature_data[,colnames(feature_data)%in%feature_selected])
  if(category == "gene"& feature_max > 50){
    rownames(feature_data)<-NULL
    feature_data<-column_to_rownames(feature_data,var = "ID")
    feature_data<-log2(feature_data)

    feature_data<-rownames_to_column(feature_data,var = "ID")
    feature_selected<-feature_manipulation(data = feature_data,feature = setdiff(colnames(feature_data),"ID"))
    feature_data<-feature_data[,colnames(feature_data)%in%c("ID",feature_selected)]
    feature_data<-as.data.frame(feature_data)
  }

  pf<-merge(pdata_group,feature_data,by = "ID",all = F)
  scale_begin<-length(colnames(pdata_group))+1
  pf[,scale_begin:ncol(pf)]<-scale( pf[,scale_begin:ncol(pf)],center = T,scale = T)


  if(!class(signature_group)=="list") stop(">>> Input must be a list.")
  ####################################################################
  for (x in 1:length(panel)) {

    index_i<- which(names(group_list)==panel[x])[1]
    group_name<-names(group_list)[index_i]

    print(paste0(">>>  Processing signature: ", group_name))

    features<-group_list[[index_i]]
    features<-features[features%in%colnames(pf)]
    if(length(features)<= 2) next

    #####################################
    #if counts of features is larger than feature limit, choose most significant variables
    if(length(features) > feature_limit){

      if(is_target_continuous==FALSE){

        eset<-pf[,colnames(pf)%in%c(group,features)]

        if(group == "group3") eset<-eset[!eset$group3=="Middle",]

        res<- batch_wilcoxon(data = eset,target = group,feature = setdiff(colnames(eset),group))
        good_features<-high_var_fea(result = res,target = "sig_names",
                                            name_padj = "p.adj",padj_cutoff = padj_cutoff,
                                            name_logfc = "statistic",logfc_cutoff  = 0,
                                            n = feature_limit/2)
      }else if(is_target_continuous==TRUE){
        eset<-pf[,colnames(pf)%in%c(target,features)]
        res<- batch_cor(data = eset,target = target,feature = setdiff(colnames(eset),target),method = "spearman")
        good_features<-high_var_fea(result = res,
                                    target = "sig_names",
                                    name_padj = "p.adj",
                                    padj_cutoff = padj_cutoff,
                                    name_logfc = "statistic",
                                    logfc_cutoff  = 0,
                                    n = feature_limit/2)
      }
      if(length(good_features)<=2){
        good_features<-res$sig_names[1:6]
        message(paste0("In panel ",group_name," : No feature is statistically significant."))
      }

      features<-features[features%in%good_features]

    }

    if(!is.null(target)){
      pf_inter<-as_tibble(pf[,c(setdiff(colnames(pdata_group),target),target,features)])

      if(target%in%colnames(pf_inter)) index_i_target<-which(colnames(pf_inter)==target)

      pf_long <- tidyr::pivot_longer(pf_inter, c(index_i_target + 1):ncol(pf_inter),
                                     names_to = "variables",values_to = "value")
    }else{
      pf_inter<-as_tibble(pf[,c("ID",group,features)])
      pf_long <- tidyr::pivot_longer(pf_inter, 3:ncol(pf_inter),
                                     names_to = "variables",values_to = "value")
    }

    # pf_long[-grep(pf_long$variables,pattern = target),]

    patterns<-tolower(gsub(patterns_to_na,pattern = "\\_",replacement = ""))
    if(tolower(group_name)%in%patterns){
      pf_long<-remove_names(input_df = pf_long,
                            variable = "variables",
                            patterns_to_na = patterns_to_na,
                            patterns_space = c("\\_"))
    }

    # pf_long$variables<-gsub(pf_long$variables,pattern = "\\_",replacement = " ")

    pf_long$variables<-substring(pf_long$variables,1,character_limit)

    if(group == "group3"){
      target_binary<-"group3"
      pf_long_group<-pf_long[!pf_long$group3=="Middle",]
    }else if(group =="best_cutoff"&!is.null(target)) {
      target_binary<-paste0(target,"_binary")
      pf_long_group<-pf_long
    }else if(group=="group2"){
      target_binary<-"group2"
      pf_long_group<-pf_long
    }else{
      target_binary<-group
      pf_long_group<-pf_long
    }

    axis_text_size<- 18 - max(nchar(pf_long$variables))/5

    color_box<-palettes(category = "box",palette = palette_box, show_col = show_col)
    ######################################################
    p <-ggboxplot(pf_long_group, x = "variables", y = "value",fill =  target_binary)+
      scale_fill_manual(values= color_box)+
      ylab(paste0(title.y))+
      # xlab("")+
      ggtitle(paste0(group_name))+
      theme_light()+
      theme(plot.title=element_text(size=rel(2),hjust=0.5),
            axis.title.y=element_text(size=rel(1.5)),
            axis.title.x = element_blank(),
            # axis.text=element_text(size=rel(2.5)),
            axis.text.x= element_text(face="plain",size=axis_text_size,
                                      angle=60,hjust = 1,color="black"),#family="Times New Roman"
            axis.text.y= element_text(face="plain",size=15,angle=0,hjust = 1,color="black"),
            axis.line=element_line(color="black",size=.5))+theme(
              legend.key.size=unit(.3,"inches"),
              legend.title = element_blank(),
              legend.position="bottom",
              legend.direction="horizontal",
              legend.justification=c(.5,.5),
              legend.box="horizontal",
              legend.box.just="top",
              legend.text=element_text(colour="black",size=10,face = "plain"))
    #################################################
    max_variables <- max(pf_long_group$value)
    group_box<-sym(target_binary)
    pp1<-p+stat_compare_means(aes(group = !!group_box,label = paste0("p = ", ..p.format..)),
                              size= 2.6, label.y = max_variables-0.3)
    pp2<-p+stat_compare_means(aes(group = !!group_box ),label = "p.signif",
                              size=7, label.y = max_variables-0.5 )
    if(show_plot&length(features)<13){
      print(pp1)
    }else if(show_plot&length(features)>13){
      print(pp2)
    }
    ###################################################
    plot_width<-length(features)*0.5+3
    plot_height<- 4 + max(nchar(pf_long_group$variables))*0.1
    ###################################################

    if(!is.null(target)){
      ggsave(pp1,filename =paste0("1-",x,"-1-",ProjectID,"-",target,"-",group_name,"-pvalue-box.pdf"),
             width = plot_width,height = plot_height,path = file_store)
      ggsave(pp2,filename =paste0("1-",x,"-2-",ProjectID,"-",target,"-",group_name,"-box.pdf"),
             width = plot_width,height = plot_height,path = file_store)
    }else{
      ggsave(pp1,filename =paste0("1-",x,"-1-",ProjectID,"-",group,"-",group_name,"-pvalue-box.pdf"),
             width = plot_width,height = plot_height,path = file_store)
      ggsave(pp2,filename =paste0("1-",x,"-2-",ProjectID,"-",group,"-",group_name,"-box.pdf"),
             width = plot_width,height = plot_height,path = file_store)
    }
    ####################################################

    #---heatmap------------------------

    colnames(pf_long_group)[which(colnames(pf_long_group)== group)]<-"target_group"
    ###################################################
    pf_long_group$value[pf_long_group$value > 2.5] = 2.5
    pf_long_group$value[pf_long_group$value < -2.5] = -2.5

    height_heatmap<-length(features)*0.15 + 3
    ####################################################
    heatmap_col<-palettes(category = "tidyheatmap",palette = palette_heatmap,show_col = show_col)
    ####################################################
    pp<-pf_long_group %>%
      group_by(target_group) %>%
      tidyHeatmap:: heatmap(
        .column = ID,
        .row = variables,
        .value = value,
        # column_title = group_name,
        # annotation = group2,
        palette_value = heatmap_col,
        show_column_names = show_heatmap_col_name)

    if(show_plot) print(pp)

     pp %>%  tidyHeatmap::save_pdf(paste0(abspath, "1-",x,"-3-",ProjectID,"-",group,"-",group_name, "-tidyheatmap.pdf"),
                                 width = 8,height = height_heatmap)
    ####################################################
    if(is_target_continuous ==TRUE& length(group_list[[index_i]])<= 20){

      #------corrplot---------------------------------
      pf_cor<-pf[,colnames(pf)%in%c(target,features)]

      if(tolower(group_name)%in%patterns){
        patterns<-tolower(gsub(patterns_to_na,pattern = "\\_",replacement = ""))
        pf_cor<-remove_names(input_df = pf_cor,
                              variable = "colnames",
                              patterns_to_na = patterns_to_na,
                              patterns_space = NULL)
      }

      bbcor <-Hmisc:: rcorr(as.matrix(pf_cor),type = "spearman")

      ###################################
      col<- palettes(category = "heatmap3",palette = palette_corplot, show_col = show_col)

      width_heatmap<-length(group_list[[index_i]])*0.75+5
      height_heatmap<-length(group_list[[index_i]])*0.75+4
      #####################################
      pdf(file  = paste0(abspath, "1-",x,"-4-",ProjectID,"-",group_name,
                         "-associated-",category, "-corplot.pdf"),
          width = width_heatmap,height = height_heatmap)
      corrplot::corrplot(bbcor$r, type="lower", order="hclust",
                         p.mat = bbcor$P, sig.level = 0.05,tl.srt=45,
                         tl.col = "black",tl.cex = 1.3,
                         addrect=2,rect.col = "black",
                         rect.lwd = 3,
                         col = colorRampPalette(col)(50))
      dev.off()
      corrplot::corrplot(bbcor$r, type="lower", order="hclust",
                         p.mat = bbcor$P, sig.level = 0.05,tl.srt=45,
                         tl.col = "black",tl.cex = 1,
                         addrect=2,rect.col = "black",
                         rect.lwd = 3,
                         col = colorRampPalette(col)(50))
      ########################################
      lab_size<- 13 - max(nchar(pf_long_group$variables))/4 #size of coefficient
      tl_cex<- 20 - max(nchar(pf_long_group$variables))/9  #size of signature name
      p<-ggcorrplot::ggcorrplot(bbcor$r, hc.order = TRUE, type = "lower", p.mat = bbcor$P, lab = TRUE,
                                pch.cex = 4.3,
                                lab_size = lab_size,
                                tl.cex =tl_cex,
                                title = names(group_list)[index_i],
                                ggtheme = ggplot2::theme_bw,
                                colors = col)+
        theme(plot.title=element_text(size=rel(2.5),hjust=0.5))
      ######################################
      ggsave(p,filename = paste0("1-",x,"-5-",ProjectID,"-",group_name,
                                 "-associated-",category, "-corplot.pdf"),
             width = 12,height = 12.8,
             path = file_store)

      ####################################
    }


  }

  features<-unique(unlist(group_list))

  if(is_target_continuous==FALSE& nlevels(pf[,group]==2)){
    eset<-pf[,colnames(pf)%in%c(group,setdiff(features,c("ID",group)))]
    res<-  batch_wilcoxon(data = eset,target = group,feature = feature_selected)
    return(res)
  }
  if(is_target_continuous==TRUE){
    eset<-pf[,colnames(pf)%in%c(target,setdiff(features,c("ID",target)))]
    res<-  batch_cor(data = eset,target = target,feature = setdiff(colnames(eset),target),method = "spearman")
    return(res)
  }

}



