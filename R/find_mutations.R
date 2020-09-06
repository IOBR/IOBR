






#' Title
#'
#' @param mutation_matrix
#' @param signature_matrix
#' @param id_signature_matrix
#' @param signature
#' @param mutation_freq
#' @param plot
#' @param method
#'
#' @return
#' @export
#'
#' @examples
find_mutations<-function(mutation_matrix, signature_matrix, id_signature_matrix = "ID", signature,mutation_freq = 0.05,plot = TRUE, method = "multi",save_path = NULL,palette = "paired3"){



  if(is.null(save_path)){
    file_name<-paste0(signature,"-relevant-mutations")
  }else{
    file_name<-save_path
  }

  if ( ! file.exists(file_name) )
    dir.create(file_name)
  abspath<-paste0(getwd(),"/",file_name,"/" )
  #######################################################

  if(max(mutation_matrix)>4){
    mutation_matrix[mutation_matrix>=3&mutation_matrix<=5]<-3
    mutation_matrix[mutation_matrix>5]<-4
  }

  #filter genes
  mut2<-mutation_matrix
  mut2[mut2>=1]<-1
  mutfreq<-data.frame(head(sort(colSums(mut2),decreasing = T),500))
  colnames(mutfreq)<-"Freq"
  index<-which(mutfreq$Freq>=dim(mut2)[1]*mutation_freq)
  index<-max(index)
  input_genes<-names(head(sort(colSums(mut2),decreasing = T),index))
  input_genes<-unique(input_genes)
  input_genes<-input_genes[!is.na(input_genes)]

  mutation_matrix<-mutation_matrix[,colnames(mutation_matrix)%in%input_genes]
  ########################################################

  colnames(signature_matrix)[which(colnames(signature_matrix)==id_signature_matrix)]<-"ID"

  mutation_matrix<-as.data.frame(mutation_matrix)
  mutation_matrix<-tibble:: rownames_to_column(mutation_matrix,var = "ID")


  sig_mut<-merge(signature_matrix,mutation_matrix,by = "ID",all = F)
  sig_mut<-tibble:: column_to_rownames(sig_mut,var = "ID")


  mytheme<-theme_light()+   ###theme_bw()  就是EXCEL的格式但是还带有表格
    theme(
      plot.title=element_text(size=rel(2.8),hjust=0.5),
      axis.title.y=element_text(size=rel(2.5)),
      axis.title.x= element_blank(),
      axis.text.x= element_text(face="plain",size=30,angle=0,color="black"),#family="Times New Roman"
      axis.text.y = element_text(face="plain",size=20,angle=90,color="black"),#family="Times New Roman"
      # panel.grid.major=element_line(color="white"), ####变成白色就是看不到了
      # panel.grid.minor=element_line(color="white"),  ####变成白色就是看不到了
      # panel.border=element_rect(color="white"),   ####变成白色就是看不到了
      axis.line=element_line(color="black",size=0.70))+theme(
        legend.key.size=unit(.3,"inches"),#图例分类符号的大小
        legend.title=element_blank(),
        legend.position="none",#"none","left","right","bottom","top",or #c(0.5,1)
        legend.direction="horizontal",# 图例排列方向#,"vertical"
        legend.justification=c(.5,.5),#"center" or two-element numeric vector
        legend.box="vertical",#"horizontal",对图例的排列方式
        legend.box.just="top",#多图例的居中方式
        legend.text=element_text(colour="black",size=10,face = "plain")#("plain", "italic", "bold", "bold.italic")
      )


  if(method == "multi"){


    input<-sig_mut[,c(signature,input_genes)]

    ##############################
    aa<-lapply(input[,input_genes], function(x) PMCMRplus::cuzickTest(input[,1]~x))
    res1<-data.frame(p.value = sapply(aa, getElement, name = "p.value"),
                    names=input_genes,
                    statistic = sapply(aa, getElement, name = "statistic"))
    res1$adjust_pvalue<-p.adjust(res1$p.value,method = "BH",n=length(res1$p.value))
    message(">>>> cuzickTest")
    print(res1[res1$p.value<0.05,])
    res1<-res1[order(res1$p.value,decreasing = F),]

    write.csv(res1,paste0(abspath,"1-cuzickTest-test-relevant-mutations.csv"))


    if(plot){

      ####################################
      #' choose top 10 genes
      top10_genes<-res1$names[1:10]
      top10_genes<- top10_genes[!is.na(top10_genes)]
      top10_genes<-as.character(top10_genes[top10_genes%in%colnames(input)])
      input_long<-input[,c(signature,top10_genes)]
      input_long<-reshape2:: melt(input_long,id.vars = 1,
                          variable.name = "Gene",
                          value.name = "mutation")
      input_long[,signature]<-as.numeric(input_long[,signature])
      input_long$mutation<-as.factor(input_long$mutation)
      ####################################

      pl<-list()
      for(i in 1:length(top10_genes)){
        gene<- top10_genes[i]
        dd<-input_long[input_long$Gene==gene,]
        pl[[i]]<-ggplot(dd, aes(x=mutation, y = !!sym(signature), fill=mutation)) +
          geom_boxplot(outlier.shape = NA,outlier.size = NA)+
          geom_jitter(width = 0.25,size= 5.9,alpha=0.75,color ="black")+
          scale_fill_manual(values=mydb:: my_palette(category = "box",palette = palette))+
          mytheme+theme(legend.position="none")+
          ggtitle(paste0(top10_genes[i]))+
          mytheme+stat_compare_means(comparisons = combn(as.character(unique(dd[,"mutation"])), 2, simplify=F),size=6)+
          stat_compare_means(size=6)
        ggsave(pl[[i]],filename = paste0(4+i,"-1-",gene,"-continue.pdf"),
               width = 4.2,height = 6.5,path = file_name)
      }
      # paste("pl","[[",1:10,"]]",sep = "",collapse=",")
      com_plot<- cowplot::plot_grid(pl[[1]],pl[[2]],pl[[3]],pl[[4]],pl[[5]],pl[[6]],pl[[7]],pl[[8]],pl[[9]],pl[[10]],
                      labels = "AUTO",ncol = 5,nrow = 2,label_size = 36)
      #####################################
      ggsave(com_plot,filename = "3-Relevant_mutations_Continue.pdf",width = 25,height = 17,path = file_name)
      #####################################
    }

    ################################
    sig_mut2<- rownames_to_column(sig_mut,var = "ID")


    patr1<-sig_mut2[,c("ID",signature)]
    part2<-sig_mut2[,input_genes]
    part2[part2>=1]<-1
    sig_mut2<-cbind(patr1,part2)
    sig_mut2<-column_to_rownames(sig_mut2,var = "ID")

    input2<-sig_mut2
    #' willcoxon test
    ##############################
    aa<-lapply(input2[,input_genes], function(x) wilcox.test(input2[,1]~x))

    res2<-data.frame(p.value = sapply(aa, getElement, name = "p.value"),
                    names=input_genes,
                    statistic = sapply(aa, getElement, name = "statistic"))

    res2$adjust_pvalue<-p.adjust(res2$p.value,method = "BH",n=length(res2$p.value))
    res2[res2$p.value<0.05,]
    res2<-res2[order(res2$p.value,decreasing = F),]


    message(">>>> wilcoxon test")
    print(res2[res2$p.value<0.05,])

    write.csv(res2, paste0(abspath,"2-Wilcoxon-test-relevant-mutations.csv"))

    result<-list("cuzick_test" = res1, "wilcoxon_test" = res2,
                 'sig_mut_data1' = input, "sig_mut_data2" = input2)


    if(plot){

      ####################################
      #' choose top 10 genes
      top10_genes<-res2$names[1:10]
      top10_genes<- top10_genes[!is.na(top10_genes)]
      top10_genes<-as.character(top10_genes[top10_genes%in%colnames(input2)])

      input_long<-input2[,c(signature,top10_genes)]
      input_long<-reshape2:: melt(input_long,id.vars = 1,
                          variable.name = "Gene",
                          value.name = "mutation")
      input_long[,signature]<-as.numeric(input_long[,signature])
      input_long$mutation<-as.factor(input_long$mutation)


      input_long$mutation<-ifelse(input_long$mutation==0,"WT","Mutated")
      ####################################

      pl<-list()
      for(i in 1:length(top10_genes)){
        gene<- top10_genes[i]
        dd<-input_long[input_long$Gene==gene,]
        pl[[i]]<-ggplot(dd, aes(x=mutation, y = !!sym(signature), fill=mutation)) +
          geom_boxplot(outlier.shape = NA,outlier.size = NA)+
          geom_jitter(width = 0.25,size=5.5,alpha=0.75,color ="black")+
          scale_fill_manual(values= mydb:: my_palette(category = "box",palette = palette))+
          mytheme+
          theme(legend.position="none")+
          ggtitle(paste0(top10_genes[i]))+
          mytheme+
          stat_compare_means(comparisons = combn(as.character(unique(dd[,"mutation"])), 2, simplify=F),size=6)
        ggsave(pl[[i]],filename = paste0(4+i,"-2-",gene,"-binary.pdf"),
               width = 4.2,height = 6.5,path = file_name)
      }
      # paste("pl","[[",1:10,"]]",sep = "",collapse=",")
      com_plot<-cowplot:: plot_grid(pl[[1]],pl[[2]],pl[[3]],pl[[4]],pl[[5]],pl[[6]],pl[[7]],pl[[8]],pl[[9]],pl[[10]],
                          labels = "AUTO",ncol = 5,nrow = 2,label_size = 36)
      #####################################
      ggsave(com_plot,filename = "4-Relevant_mutations_binary.pdf",width = 22,height = 17,path = file_name)
      #####################################
    }


  }else{
    ################################
    patr1<-sig_mut[,signature]
    part2<-sig_mut[,input_genes]
    part2[part2>=1]<-1
    sig_ms2<-cbind(patr1,part2)

    input2<-sig_ms2
    #' willcoxon test
    ##############################
    aa<-lapply(input2[,input_genes], function(x) wilcox.test(input2[,1]~x))

    res<-data.frame(p.value = sapply(aa, getElement, name = "p.value"),
                     names=input_genes,
                     statistic = sapply(aa, getElement, name = "statistic"))

    res$adjust_pvalue<-p.adjust(res$p.value,method = "BH",n=length(res$p.value))
    res[res$p.value<0.05,]
    res<-res[order(res$p.value,decreasing = F),]

    write.csv(res, paste0(abspath,"1-Wilcoxon-test-relevant-mutations.csv"))


    result<-list("wilcoxon_test" = res,"sig_mut_data" = input2)
    #################################

    if(plot){

      ####################################
      #' choose top 10 genes
      top10_genes<-res$names[1:10]
      top10_genes<-as.character(top10_genes[top10_genes%in%colnames(input2)])

      input_long<-input2[,c(signature,top10_genes)]
      input_long<-reshape2:: melt(input_long,id.vars = 1,
                          variable.name = "Gene",
                          value.name = "mutation")
      input_long[,signature]<-as.numeric(input_long[,signature])
      input_long$mutation<-as.factor(input_long$mutation)

      input_long$mutation<-ifelse(input_long$mutation==0,"WT","Mutated")
      ####################################
      pl<-list()
      for(i in 1:length(top10_genes)){
        gene<- top10_genes[i]
        dd<-input_long[input_long$Gene==gene,]
        pl[[i]]<-ggplot(dd, aes(x=mutation, y = !!sym(signature), fill=mutation)) +
          geom_boxplot(outlier.shape = NA,outlier.size = NA)+
          geom_jitter(width = 0.25,size=5.5,alpha=0.75,color ="black")+
          scale_fill_manual(values= mydb::my_palette(category = "box",palette = palette))+
          mytheme+theme(legend.position="none")+
          ggtitle(paste0(top10_genes[i]))+
          mytheme+
          stat_compare_means(comparisons = combn(as.character(unique(dd[,"mutation"])), 2, simplify=F),size=6)

        ggsave(pl[[i]],filename = paste0(2+i,"-1-",gene,"-binary.pdf"),
               width = 4.2,height = 6.5,path = file_name)
      }
      # paste("pl","[[",1:10,"]]",sep = "",collapse=",")
      com_plot<-cowplot:: plot_grid(pl[[1]],pl[[2]],pl[[3]],pl[[4]],pl[[5]],pl[[6]],pl[[7]],pl[[8]],pl[[9]],pl[[10]],
                          labels = "AUTO",ncol = 5,nrow = 2,label_size = 36)
      #####################################
      ggsave(com_plot,filename = "2-Relevant_mutations_binary.pdf",width = 22,height = 17,path = file_name)
      #####################################
    }


  }

  return(result)


}
