


#' List of supported signature score calculation methods
#'
#' The methods currently supported are
#' `PCA`, `ssGSEA`, `z-score`, `Integration`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export

signature_score_calculation_methods= c("PCA" = "pca",
                                      "ssGSEA" = "ssgsea",
                                      "z-score" = "zscore",
                                      "Integration" = "integration")
##########################################




#' Calculating signature score using PCA method
#'
#' @param pdata phenotype data of input sample;
#' if phenotype data is NULL, create a data frame with `Index` and `ID` contain column names of eset
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures;
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set
#' @param column_of_sample  Defines in which column of pdata the sample identifier can be found.
#'
#' @author Dongqiang Zeng
#' @return matrix;pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#'
#' @examples
calculate_sig_score_pca<-function(pdata,
                                  eset,
                                  signature,
                                  mini_gene_count,
                                  column_of_sample){

  message(paste0("\n", ">>> Calculating signature score with PCA method"))

  #creat pdata
  if(is.null(pdata)){
    pdata<-data.frame("Index" = 1:length(colnames(eset)),"ID" = colnames(eset))
  }

  #match phenotype data and gene expression set
  ###########################
  colnames(pdata)[which(colnames(pdata)==column_of_sample)]<-"ID"
  pdata<-pdata[pdata$ID%in%colnames(eset),]
  eset<-eset[,colnames(eset)%in%pdata$ID]
  eset<-eset[,match(pdata$ID,colnames(eset))]
  ###########################

  #normalization
  if(max(eset)>100) eset<-log2(eset+1)
  eset<-scale(eset,center = T,scale = T)
  ###########################

  #filter signatures
  if(mini_gene_count<=2) mini_gene_count <- 2
  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= mini_gene_count]
  ###########################

  #calculating signature score
  goi <- names(signature)
  ###########################
  for (sig in goi) {
    pdata[, sig] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop=FALSE]
    pdata[, sig] <- sigScore(tmp,methods = "PCA")
  }
  if ("TMEscoreA_CIR"%in%goi &"TMEscoreB_CIR" %in% goi) {
    pdata[,"TMEscore_CIR"]<-pdata[,"TMEscoreA_CIR"]-pdata[,"TMEscoreB_CIR"]
  }
  if ("TMEscoreA_plus"%in%goi &"TMEscoreB_plus" %in% goi) {
    pdata[,"TMEscore_plus"]<-pdata[,"TMEscoreA_plus"]-pdata[,"TMEscoreB_plus"]
  }
  return(pdata)
}
###################################################






#' Calculating signature score using z-score method
#'
#' @param pdata phenotype data of input sample;
#' if phenotype data is NULL, create a data frame with `Index` and `ID` contain column names of eset
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set;
#' @param column_of_sample  Defines in which column of pdata the sample identifier can be found
#'
#' @author Dongqiang Zeng
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#'
#' @examples
#'
calculate_sig_score_zscore<-function(pdata, eset, signature,
                                     mini_gene_count,
                                     column_of_sample){

  message(paste0("\n", ">>> Calculating signature score with z-score method"))

  #creat pdata
  if(is.null(pdata)){
    pdata<-data.frame("Index" = 1:length(colnames(eset)),"ID" = colnames(eset))
  }
  #match phenotype data and gene expression set
  ###########################
  colnames(pdata)[which(colnames(pdata)==column_of_sample)]<-"ID"
  pdata<-pdata[pdata$ID%in%colnames(eset),]
  eset<-eset[,colnames(eset)%in%pdata$ID]
  eset<-eset[,match(pdata$ID,colnames(eset))]
  ###########################
  #normalization
  if(max(eset)>100) eset<-log2(eset+1)
  eset<-scale(eset,center = T,scale = T)
  ###########################
  ###########################
  if(mini_gene_count<=2) mini_gene_count <- 2
  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= mini_gene_count]
  ###########################
  #calculating signature score
  goi <- names(signature)
  for (sig in goi) {
    pdata[, sig] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop=FALSE]
    pdata[, sig] <- sigScore(tmp,methods = "zscore")
  }
  if ("TMEscoreA_CIR"%in%goi &"TMEscoreB_CIR" %in% goi) {
    pdata[,"TMEscore_CIR"]<-pdata[,"TMEscoreA_CIR"]-pdata[,"TMEscoreB_CIR"]
  }
  if ("TMEscoreA_plus"%in%goi &"TMEscoreB_plus" %in% goi) {
    pdata[,"TMEscore_plus"]<-pdata[,"TMEscoreA_plus"]-pdata[,"TMEscoreB_plus"]
  }
  return(pdata)
}
###################################################




#' Calculating signature score using ssGSEA method
#'
#' @param pdata phenotype data of input sample;
#' if phenotype data is NULL, create a data frame with `Index` and `ID` contain column names of eset
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set; default is 5;
#' the minimal gene count for ssGSEA methods should larger than 5 for the robustness of the calculation
#' @param column_of_sample  Defines in which column of pdata the sample identifier can be found.
#'
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @import GSVA
#' @import tibble
#' @author Dongqiang Zeng
#' @examples
calculate_sig_score_ssgsea<-function(pdata, eset, signature,
                                     mini_gene_count,
                                     column_of_sample){

  message(paste0("\n", ">>> Calculating signature score with ssGSEA method"))

  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= mini_gene_count]
  ###########################
  if(mini_gene_count<=5) mini_gene_count <- 5
  ############################

  #creat pdata
  if(is.null(pdata)){
    pdata<-data.frame("Index" = 1:length(colnames(eset)),"ID" = colnames(eset))
  }
  ##############################
  res <- GSVA:: gsva(as.matrix(eset),
                     signature,
                     method="ssgsea",
                     kcdf="Gaussian",
                     min.sz= mini_gene_count,
                     ssgsea.norm=T)
  res<-as.data.frame(t(res))
  res<-rownames_to_column(res,var = "ID")

  if ("TMEscoreA_CIR"%in%colnames(res) & "TMEscoreB_CIR"%in%colnames(res)) {
    res[,"TMEscore_CIR"]<-res[,"TMEscoreA_CIR"] - res[,"TMEscoreB_CIR"]
  }

  if ("TMEscoreA_plus"%in%colnames(res) & "TMEscoreB_plus"%in%colnames(res)) {
    res[,"TMEscore_plus"]<-res[,"TMEscoreA_plus"] - res[,"TMEscoreB_plus"]
  }

  pdata<-merge(pdata,res,by ="ID",all.x = T,all.y = F)
  return(pdata)
}
###################################################






#' Calculating signature score using Integration method
#'
#' @param pdata phenotype data of input sample;
#' if phenotype data is NULL, create a data frame with `Index` and `ID` contain column names of eset
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set; default is 5;
#' the minimal gene count for ssGSEA methods should larger than 5 for the robustness of the calculation
#' @param column_of_sample  Defines in which column of pdata the sample identifier can be found.
#'
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @import GSVA
#' @import tibble
#' @author Dongqiang Zeng
#' @examples
#'
calculate_sig_score_integration<-function(pdata, eset, signature,
                                          mini_gene_count = 2,
                                          column_of_sample){
  message(paste0("\n", ">>> Calculating signature score with Integration methods"))

  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= mini_gene_count]
  ###########################
  #creat pdata if NULL
  if(is.null(pdata)){
    pdata<-data.frame("Index" = 1:length(colnames(eset)),"ID" = colnames(eset))
  }
  #match phenotype data and gene expression set
  ###########################
  colnames(pdata)[which(colnames(pdata)==column_of_sample)]<-"ID"
  pdata<-pdata[pdata$ID%in%colnames(eset),]
  eset<-eset[,colnames(eset)%in%pdata$ID]
  eset<-eset[,match(pdata$ID,colnames(eset))]
  ###########################
  #normalization
  if(max(eset)>100) {
    eset1<-log2(eset+1)
  }else{eset1<-eset}
  ##########################
  eset1<-scale(eset1,center = T,scale = T)
  message(paste0("\n", ">>>Step 1: Calculating signature score using PCA method"))
  goi <- names(signature)
  ###########################
  for (sig in goi) {
    pdata[, paste(sig,"_PCA",sep = "")] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset1)]
    tmp <- eset1[genes, , drop=FALSE]
    pdata[, paste(sig,"_PCA",sep = "")] <- sigScore(tmp,methods = "PCA")
  }
  ############################
  if ("TMEscoreA_CIR"%in%goi &"TMEscoreB_CIR" %in% goi) {
    pdata[,"TMEscore_CIR_PCA"]<-pdata[,"TMEscoreA_CIR_PCA"]-pdata[,"TMEscoreB_CIR_PCA"]
  }
  if ("TMEscoreA_plus"%in%goi &"TMEscoreB_plus" %in% goi) {
    pdata[,"TMEscore_plus_PCA"]<-pdata[,"TMEscoreA_plus_PCA"]-pdata[,"TMEscoreB_plus_PCA"]
  }
  ############################
  message(paste0("\n", ">>>Step 2: Calculating signature score using z-score method"))
  for (sig in goi) {
    pdata[, paste(sig,"_zscore",sep = "")] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset1)]
    tmp <- eset1[genes, , drop=FALSE]
    pdata[, paste(sig,"_zscore",sep = "")] <- sigScore(tmp,methods = "z-score")
  }
  if ("TMEscoreA_CIR"%in%goi &"TMEscoreB_CIR" %in% goi) {
    pdata[,"TMEscore_CIR_zscore"]<-pdata[,"TMEscoreA_CIR_zscore"]-pdata[,"TMEscoreB_CIR_zscore"]
  }
  if ("TMEscoreA_plus"%in%goi &"TMEscoreB_plus" %in% goi) {
    pdata[,"TMEscore_plus_zscore"]<-pdata[,"TMEscoreA_plus_zscore"]-pdata[,"TMEscoreB_plus_zscore"]
  }
  ############################
  message(paste0("\n", ">>>Step 3: Calculating signature score using ssGSEA method"))
  res <- gsva(as.matrix(eset), signature,
                  method="ssgsea",
                  kcdf="Gaussian", #标准化好的所有基因表达矩阵都推荐使用Gaussian
                  min.sz=5,
                  ssgsea.norm=T)
  res<-as.data.frame(t(res))
  if ("TMEscoreA_CIR"%in% colnames(res) & "TMEscoreB_CIR"%in%colnames(res)) {
    res[,"TMEscore_CIR"]<-res[,"TMEscoreA_CIR"] - res[,"TMEscoreB_CIR"]
  }
  if ("TMEscoreA_plus"%in% colnames(res) & "TMEscoreB_plus"%in%colnames(res)) {
    res[,"TMEscore_plus"]<-res[,"TMEscoreA_plus"] - res[,"TMEscoreB_plus"]
  }
  #############################
  colnames(res)<-paste0(colnames(res),"_ssGSEA")
  res<-rownames_to_column(res,var = "ID")
  pdata<-merge(pdata,res,by="ID",all.x = T,all.y = T)
  return(pdata)
}






#' Calculating signature score  on a gene expression dataset
#' @param pdata phenotype data of input sample
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures;
#' such as `signature_tme`, `signature_metabolism`,`signature_collection`, `go_bp`,`kegg`,`hallmark`
#' @param method he methods currently supported are `pca`, `ssgsea`, `zscore`,`integration`
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set;
#' default is 2 for PCA and z score funciion
#' @param column_of_sample Defines in which column of pdata the sample identifier can be found.
#' @param print_gene_propotion logical, print the propotion of signature genes in gene matrix
#' @param print_filtered_signatures logical, print filtered signatures has gene count less than minical gene count
#' @param ...
#'
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @author Dongqiang Zeng
#' @examples
#'
calculate_sig_score<-function(pdata = NULL, eset,
                              signature = signature_collection,
                              method = "pca",
                              mini_gene_count = 3,
                              column_of_sample = "ID",
                              print_gene_propotion = FALSE,
                              print_filtered_signatures = FALSE,...){
  ########################################
  if (print_gene_propotion) {
    message(lapply(signature,function(x) summary(x%in%rownames(eset))))
  }
  ########################################
  if (print_filtered_signatures) {
    filtered_signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))<= mini_gene_count]
    message(paste0("\n"," The number of filtered signatures is: ", length(filtered_signature)))
    if (length(filtered_signature)>=10) {
      message(paste("\n", "10 of signatures with gene count less than ", mini_gene_count, " : "))
      message(paste( " ", names(filtered_signature)[1:10],collapse = "\n"))
      message(paste0("\n", "You can use this function to get all filtered sigantures:", "\n",
                     "signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))<= mini_gene_count]"))
    }else if(length(filtered_signature)>0 & length(filtered_signature) < 10){
      message(paste( names(filtered_signature)[1:length(filtered_signature)],collapse = "\n"))
    }else if(length(filtered_signature)<=0){
      message(paste0( "\n", "All signatures are eligible for calculating signature score"))
    }
  }
  ##########################################
  # run selected method
  res = switch(method,
               pca = calculate_sig_score_pca(pdata, eset,
                                             signature = signature,
                                             mini_gene_count = mini_gene_count,
                                             column_of_sample = column_of_sample,...),
               ssgsea = calculate_sig_score_ssgsea(pdata, eset,
                                                   signature = signature,
                                                   mini_gene_count = mini_gene_count,
                                                   column_of_sample = column_of_sample,...),
               zscore = calculate_sig_score_zscore(pdata,eset,
                                                   signature = signature,
                                                   mini_gene_count = mini_gene_count,
                                                   column_of_sample = column_of_sample,...),
               integration = calculate_sig_score_integration(pdata,eset,
                                                        signature = signature,
                                                        mini_gene_count = mini_gene_count,
                                                        column_of_sample = column_of_sample,...))
  return(res)
}








