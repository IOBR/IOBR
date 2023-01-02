





#' Make mutation matrix using maf data
#'
#' @param maf file name of maf data
#' @param mut_data if maf data had beed loaded, this argument will be applied
#' @param isTCGA logical variable,
#' @param category if category is `multi`, four mutation matrix will be generated: all, snp, indel and frameshift
#' @param Tumor_Sample_Barcode column name of tumor sample
#' @param Hugo_Symbol column name of gene symbol
#' @param Variant_Classification column name of variant classification: Frame_Shift_Del
#' @param Variant_Type column name of variant type: SNP, INS, DEL
#'
#' @return
#' @export
#' @author Dongqian Zeng
#'
#' @examples
make_mut_matrix<-function(maf = NULL, mut_data = NULL, isTCGA = TRUE, category = "multi", Tumor_Sample_Barcode = "Tumor_Sample_Barcode", Hugo_Symbol = "Hugo_Symbol", Variant_Classification = "Variant_Classification", Variant_Type = "Variant_Type"){


  if(!is.null(maf)){

    mut_maf<-maftools:: read.maf(maf = maf,useAll = TRUE,isTCGA = isTCGA)

    print(summary(mut_maf@data$Variant_Classification))
    print(summary(mut_maf@data$Variant_Type))
    #########################################

    mut<-mut_maf@data
  }else{

    mut<-as.data.frame(mut_data)
    colnames(mut)[which(colnames(mut)==Tumor_Sample_Barcode)]<-"Tumor_Sample_Barcode"
    colnames(mut)[which(colnames(mut)==Hugo_Symbol)]<-"Hugo_Symbol"
    colnames(mut)[which(colnames(mut)==Variant_Classification)]<-"Variant_Classification"

    if(!Variant_Type%in%colnames(mut)){

      if("filter"%in%colnames(mut)) mut<-mut[mut$filter=="PASS",]
      mut$Variant_Type<-mut$Variant_Classification
      mut<-mut[!grepl(mut$Variant_Type,pattern = "synonymous"),]

      mut$Variant_Type[stringr::str_detect(tolower(mut$Variant_Type),pattern = "missense")]<-"SNP"
      mut$Variant_Type[stringr::str_detect(tolower(mut$Variant_Type),pattern = "frameshift")]<-"Frame_Shift"
      mut$Variant_Type[stringr::str_detect(tolower(mut$Variant_Type),pattern = "insert")]<-"INS"
      mut$Variant_Type[stringr::str_detect(tolower(mut$Variant_Type),pattern = "delet")]<-"DEL"

      mut<-mut[mut$Variant_Type%in%c("SNP","INS","DEL","Frame_Shift"),]
    }else{
      colnames(mut)[which(colnames(mut)==Variant_Type)]<-"Variant_Type"
    }

  }

  mut<-mut[,c("Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","Variant_Type")]

  if(category == "multi"){

    mut_all <-reshape2:: dcast(mut, Hugo_Symbol~Tumor_Sample_Barcode, value.var="Variant_Classification")
    mut_all<-tibble:: column_to_rownames(mut_all,var = "Hugo_Symbol")
    mut_all<-as.matrix(mut_all)
    mut_all<-matrix(as.numeric(mut_all), dim(mut_all), dimnames = dimnames(mut_all))

    mut_snp<-mut[mut$Variant_Type=="SNP",]
    mut_snp <-reshape2:: dcast(mut_snp, Hugo_Symbol~Tumor_Sample_Barcode, value.var="Variant_Classification")
    mut_snp<-tibble:: column_to_rownames(mut_snp,var = "Hugo_Symbol")
    mut_snp<-as.matrix(mut_snp)
    mut_snp<-matrix(as.numeric(mut_snp), dim(mut_snp), dimnames = dimnames(mut_snp))

    mut_indel<-mut[mut$Variant_Type=="DEL"|mut$Variant_Type=="INS",]
    mut_indel <- reshape2::dcast(mut_indel, Hugo_Symbol~Tumor_Sample_Barcode, value.var="Variant_Classification")
    mut_indel<-tibble:: column_to_rownames(mut_indel,var = "Hugo_Symbol")
    mut_indel<-as.matrix(mut_indel)
    mut_indel<-matrix(as.numeric(mut_indel), dim(mut_indel), dimnames = dimnames(mut_indel))

    mut_frameshift<-mut[grepl(tolower(mut$Variant_Classification),pattern = "frame"),]
    mut_frameshift <- reshape2::dcast(mut_frameshift, Hugo_Symbol~Tumor_Sample_Barcode, value.var="Variant_Classification")
    mut_frameshift<-tibble:: column_to_rownames(mut_frameshift,var = "Hugo_Symbol")
    mut_frameshift<-as.matrix(mut_frameshift)
    mut_frameshift<-matrix(as.numeric(mut_frameshift), dim(mut_frameshift), dimnames = dimnames(mut_frameshift))
    ####################################
    ####################################
    mut_list = list(all = t(mut_all),snp=t(mut_snp),indel=t(mut_indel),frameshift = t(mut_frameshift))
    mut_list = ComplexHeatmap:: unify_mat_list(mut_list)
  }else{

    mut<-mut[,c("Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","Variant_Type")]
    mut_all <-reshape2:: dcast(mut, Hugo_Symbol~Tumor_Sample_Barcode, value.var="Variant_Classification")
    mut_all<-tibble:: column_to_rownames(mut_all,var = "Hugo_Symbol")
    mut_all<-as.matrix(mut_all)
    mut_all<-matrix(as.numeric(mut_all), dim(mut_all), dimnames = dimnames(mut_all))
    mut_list = list(all = t(mut_all))
  }
  return(mut_list)
}
