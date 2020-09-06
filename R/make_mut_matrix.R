





#' Title
#'
#' @param maf 
#' @param maf_data 
#' @param isTCGA 
#' @param category 
#' @param Tumor_Sample_Barcode 
#' @param Hugo_Symbol 
#' @param Variant_Classification 
#' @param Variant_Type 
#'
#' @return
#' @export
#'
#' @examples
make_mut_matrix<-function(maf = NULL, maf_data = NULL, isTCGA = TRUE, category = "multi",Tumor_Sample_Barcode = "Tumor_Sample_Barcode", Hugo_Symbol = "Hugo_Symbol", Variant_Classification = "Variant_Classification", Variant_Type = "Variant_Type"){
  
  
  if(!is.null(maf)){
    
    mut_maf<-maftools:: read.maf(maf = maf,useAll = FALSE,isTCGA = isTCGA)
    
    # summary(as.factor(substring(mut_maf@variants.per.sample$Tumor_Sample_Barcode,14,16)))
    # summary(duplicated(substring(mut_maf@variants.per.sample$Tumor_Sample_Barcode,1,12)))
    #########################################
    # dim(mut_maf@data)
    print(summary(mut_maf@data$Variant_Classification))
    print(summary(mut_maf@data$Variant_Type))
    #########################################
    
    mut<-mut_maf@data
  }else{
    
    mut<-maf_data
    colnames(mut)[which(colnames(mut)==Tumor_Sample_Barcode)]<-"Tumor_Sample_Barcode"
    colnames(mut)[which(colnames(mut)==Hugo_Symbol)]<-"Hugo_Symbol"
    colnames(mut)[which(colnames(mut)==Variant_Classification)]<-"Variant_Classification"
    colnames(mut)[which(colnames(mut)==Variant_Type)]<-"Variant_Type"
    
  }
  
  mut<-mut[,c("Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","Variant_Type")]
  
  if(category == "multi"){
    
    mut_all <-reshape2:: dcast(mut, Hugo_Symbol~Tumor_Sample_Barcode, value.var="Variant_Classification")
    mut_all<-tibble:: column_to_rownames(mut_all,var = "Hugo_Symbol")
    mut_all<-as.matrix(mut_all)
    mut_all<-matrix(as.numeric(mut_all), dim(mut_all), dimnames = dimnames(mut_all))
    
    mut_snp<-mut[mut$Variant_Type=="SNP",];dim(mut_snp)
    mut_snp <-reshape2:: dcast(mut_snp, Hugo_Symbol~Tumor_Sample_Barcode, value.var="Variant_Classification")
    mut_snp<-tibble:: column_to_rownames(mut_snp,var = "Hugo_Symbol")
    mut_snp<-as.matrix(mut_snp)
    mut_snp<-matrix(as.numeric(mut_snp), dim(mut_snp), dimnames = dimnames(mut_snp))
    
    mut_indel<-mut[mut$Variant_Type=="DEL"|mut$Variant_Type=="INS",];dim(mut_indel)
    mut_indel <- reshape2::dcast(mut_indel, Hugo_Symbol~Tumor_Sample_Barcode, value.var="Variant_Classification")
    mut_indel<-tibble:: column_to_rownames(mut_indel,var = "Hugo_Symbol")
    mut_indel<-as.matrix(mut_indel)
    mut_indel<-matrix(as.numeric(mut_indel), dim(mut_indel), dimnames = dimnames(mut_indel))
    
    mut_frameshift<-mut[mut$Variant_Classification=="Frame_Shift_Del"|mut$Variant_Classification=="Frame_Shift_Ins",]
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
