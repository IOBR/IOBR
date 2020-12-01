




#' Collection of tumor microenvironment cell fraction deconvolution methods.
#'
#' The methods currently supported are
#' `mcp_counter`, `epic`, `xcell`, `cibersort`, `cibersort_abs`, `ips`, `estimate`, `svm`,`lsei`,`timer`,`quantiseq`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
tme_deconvolution_methods = c("MCPcounter"="mcpcounter",
                              "EPIC"="epic",
                              "xCell"="xcell",
                              "CIBERSORT"="cibersort",
                              "CIBERSORT Absolute"="cibersort_abs",
                              "IPS" = "ips",
                              "ESTIMATE" = "estimate",
                              "SVM" = "svm",
                              "lsei" = "lsei",
                              "TIMER" = "timer",
                              "quanTIseq" = "quantiseq")
############################################




#' Decoding immune microenvironment using xCell
#'
#' @param eset expression set with genes at row, sample ID at column
#' @param project project name used to distinguish different data sets
#' @param array transcrptomic data type
#' @return xCell with immune cell fractions
#' @export
#' @importFrom tibble rownames_to_column
#' @author Dongqiang Zeng
#' @examples
#' xcell_result<-deconvo_xcell(eset = eset_stad,project = "TCGA-STAD")
#'
deconvo_xcell<-function(eset,project = NULL,arrays){

  message(paste0("\n", ">>> Running ", "xCell"))
  #normalize gene expression matrix
  # if(max(eset)>100) eset<-log2(eset)
  # data("xCell.data")
  rnaseq = !arrays
  res<- xCellAnalysis(eset,rnaseq = rnaseq)
  res<-as.data.frame(t(res))
  ###########################################
  colnames(res)<-gsub(colnames(res),pattern = "\\ ",replacement = "\\_")
  colnames(res)<-gsub(colnames(res),pattern = "\\ ",replacement = "\\_")
  colnames(res)<-paste0(colnames(res),"_xCell")

  if(!is.null(project)){
    res$ProjectID<-project
    res<-res[,c(ncol(res),1:ncol(res)-1)]
  }
  res<-rownames_to_column(res,var = "ID")
  return(res)
}



#' Estimating immune microenvironment using MCP-counter
#'
#' @param eset expression set with genes at row, sample ID at column
#' @param project project name used to distinguish different datasets
#' @return MCP-counter with immune cell fractions
#' @export
#' @importFrom MCPcounter MCPcounter.estimate
#' @importFrom tibble rownames_to_column
#' @author Dongqiang Zeng
#' @examples
#' mcp_result<-deconvo_mcpcounter(eset = eset_stad,project = "TCGA-STAD")
#'
deconvo_mcpcounter<-function(eset,project = NULL){

  message(paste0("\n", ">>> Running ", "MCP-counter"))
  #normalize gene expression matrix
  # if(max(eset)>100) eset<-log2(eset)

  res<-MCPcounter::MCPcounter.estimate(eset,
                                       featuresType = "HUGO_symbols",
                                       probesets= mcp_probesets,
                                       genes= mcp_genes)
  res<-as.data.frame(t(res))
  ####################################
  colnames(res)<-gsub(colnames(res),pattern = "\\.",replacement = "\\_")
  colnames(res)<-gsub(colnames(res),pattern = "\\ ",replacement = "\\_")
  colnames(res)<-paste(colnames(res),"_MCPcounter",sep = "")

  if(!is.null(project)){
    res$ProjectID<-project
    res<-res[,c(ncol(res),1:ncol(res)-1)]
  }

  res<-rownames_to_column(res,var = "ID")
  return(res)
  ###################################
}





#' Estimating immune microenvironment using EPIC: FOR RNAseq mostly
#'
#' @param eset expression set with genes at row, sample ID at column
#' @param project project name used to distinguish different data sets
#' @param tumor is input sample tumor or normal
#' @return EPIC with immune cell fractions
#' @export
#' @importFrom EPIC EPIC
#' @importFrom tibble rownames_to_column
#' @author Dongqiang Zeng
#' @examples
#' epic_result<-deconvo_epic(eset = eset_stad,project = "TCGA-STAD",tumor = TRUE)
#'
deconvo_epic<-function(eset,project = NULL,tumor){

  message(paste0("\n", ">>> Running ", "EPIC"))

  # mRNA_cell = NULL
  # if(!scale_eset) mRNA_cell = c("default"=1.)
  ###################################

  ref = ifelse(tumor, "TRef", "BRef")
  ##############################
  out <- EPIC::EPIC(bulk = eset,reference = ref,mRNA_cell = NULL, scaleExprs =TRUE)
  res<-as.data.frame((out$cellFractions))
  ####################################
  colnames(res)<-gsub(colnames(res),pattern = "\\.",replacement = "\\_")
  colnames(res)<-gsub(colnames(res),pattern = "\\ ",replacement = "\\_")
  colnames(res)<-paste(colnames(res),"_EPIC",sep = "")
  res<-as.data.frame(res)

  if(!is.null(project)){
    res$ProjectID<-project
    res<-res[,c(ncol(res),1:ncol(res)-1)]
  }
  res<-rownames_to_column(res,var = "ID")
  return(res)
  ###################################
}





#' Decoding immune microenvironment using CIBERSORT
#'
#' CIBERSORT is only freely available to academic users.
#' A license an the binary can be obtained from https://cibersort.stanford.edu.
#'
#' @param eset expression set with gene symbol at rowname, sample ID at column
#' @param project project name used to distinguish different datasets
#' @param arrays logical: Runs methods in a mode optimized for microarray data.
#' @param absolute logical: Runs CIBERSORT in absolute mode
#' @param perm permutation to run CIBERSORT
#' @return cibersrot with immune cell fractions
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' cibersort_result<-deconvo_cibersort(eset = eset_stad,project = "TCGA-STAD",arrays = TRUE,absolute = FALSE, perm = 500)


deconvo_cibersort<-function(eset, project = NULL, arrays,  perm = 1000, absolute = FALSE, abs_method = "sig.score"){

  if(absolute){
    message(paste0("\n", ">>> Running ", "CIBERSORT in absolute mode"))
  }else{
    message(paste0("\n", ">>> Running ", "CIBERSORT"))
  }

  eset<-as.data.frame(eset)
  ##############################

  # the authors reccomend to disable quantile normalizeation for RNA seq.
  # (see CIBERSORT website).
  quantile_norm = arrays
  ##############################
  res<-CIBERSORT(sig_matrix = lm22,
                 mixture_file = eset,
                 perm = perm,
                 QN = quantile_norm,
                 absolute = absolute,
                 abs_method = abs_method)
  ###############################
  colnames(res)<-gsub(colnames(res),pattern = "\\.",replacement = "\\_")
  colnames(res)<-gsub(colnames(res),pattern = "\\ ",replacement = "\\_")
  colnames(res)<-paste(colnames(res),"_CIBERSORT",sep = "")
  res<-as.data.frame(res)


  if(!is.null(project)){
    res$ProjectID<-project
    res<-res[,c(ncol(res),1:ncol(res)-1)]
  }

  res<-rownames_to_column(res,var = "ID")
  return(res)
  ###################################
}




#' Calculating immune phenotype score using IPS
#'
#' @param eset expression set with genes at row, sample ID at column
#' @param project project name used to distinguish different data sets
#' @return IPS data frame
#' @export
#' @import ggplot2
#' @import cowplot
#' @import grid
#' @author Dongqiang Zeng
#' @examples
#' ips_result<-deconvo_ips(eset = eset_stad,project = "TCGA-STAD")
#'
deconvo_ips<-function(eset,project = NULL,plot){

  message(paste0("\n", ">>> Running ", "Immunophenoscore"))
  #normalize gene expression matrix
  # if(max(eset)>100) eset<-log2(eset+1)
  ##############################
  res<-IPS_calculation(project = project,eset,plot)
  ####################################
  # colnames(res)<-gsub(colnames(res),pattern = "\\.",replacement = "\\_")
  colnames(res)<-paste(colnames(res),"_IPS",sep = "");colnames(res)
  res<-as.data.frame(res)

  if(!is.null(project)){
    res$ProjectID<-project
    res<-res[,c(ncol(res),1:ncol(res)-1)]
  }

  res<-rownames_to_column(res,var = "ID")
  return(res)
}

###########################################


#' Calculation of stromal, immune, and ESTIMATE scores
#'
#' @param eset expression set with genes at row, sample ID at column
#' @param project project name used to distinguish different datasets
#' @param platform character string indicating platform type. Defaults to "affymetrix"
#' @importFrom estimate estimateScore
#' @importFrom tibble rownames_to_column
#' @author Dongqiang Zeng
#' @return
#' @export
#'
#' @examples
deconvo_estimate<-function(eset, project = NULL,platform = "affymetrix"){

  message(paste0("\n", ">>> Running ", "ESTIMATE"))
  eset<-as.data.frame(eset)
  eset<-rownames_to_column(eset,var = "symbol")
  sampleData<-paste0(project,"-eset.txt");sampleData
  write.table(eset,sampleData,sep = "\t",row.names = F,quote = F)
  ########################################
  filterCommonGenes(input.f= sampleData,
                    output.f= paste0(project,"_Tumor_purity.gct"),
                    id="GeneSymbol")
  #delete-data-after-input
  file.remove(paste0(project,"-eset.txt"))
  ################################
  estimate::estimateScore(input.ds = paste0(project,"_Tumor_purity.gct"),
                output.ds= paste0(project,"_Tumor_estimate_score.gct"),
                platform= platform)
  file.remove(paste0(project,"_Tumor_purity.gct"))
  #################################
  scores=read.table(paste0(project,"_Tumor_estimate_score.gct"),skip = 2,header = T)
  file.remove(paste0(project,"_Tumor_estimate_score.gct"))
  ################################
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  colnames(scores)<-paste0(colnames(scores),"_estimate")

  if(!is.null(project)){
    scores$ProjectID<-project
    scores<-scores[,c(ncol(scores),1:ncol(scores)-1)]
  }

  scores<-rownames_to_column(as.data.frame(scores),var = "ID")
  scores$ID<-gsub(scores$ID,pattern = "\\.",replacement = "-")

  return(scores)
}



#' Estimate the fraction of cell types using defined reference genes
#'
#' @param eset expression data with matched gene id of reference
#' @param project project name used to distinguish different datasets
#' @param arrays a logical value. If TRUE, the columns of the input data will be normalized to have the same quantiles.
#' @param method deconvolution method. must be "svm" or "lsei"
#' @param perm  permutation to run CIBERSORT
#' @param reference immune cell gene matrix; eg lm22, lm6 or can be generate using generateRef/generateRef_rnaseq
#' @param absolute.mode absolute.mode, default is FALSE
#' @param abs.method abs.method, default is sig.score
#' @param scale_reference  a logical value indicating whether the reference be scaled or not. If TRUE, the value in reference file will be centered and scaled in row direction.
#'
#' @author Dongqiang Zeng
#' @author Rongfang Shen
#' @return
#' @export
#'
#' @examples
deconvo_ref<-function(eset,project = NULL,arrays,method = "svm",perm,
                      reference, scale_reference, absolute.mode = FALSE, abs.method = "sig.score"){

  if (length(intersect(rownames(eset), rownames(reference))) == 0){
    stop("None identical gene between eset and reference had been found.
         Check your eset using: intersect(rownames(eset), rownames(reference))")
  }
  # reccomend to disable quantile normalizeation for RNA seq.
  quantile_norm = arrays
  ##############################

  if(method=="svm"){
    message(paste0("\n", ">>> Running ", "cell estimation in SVM mode"))

   # res<- deconvo_constru_sig(eset = eset,reference = reference,scale_reference = scale_reference,
   #                           method = "svm",perm = perm,arrays = arrays,
   #                           absolute.mode  = absolute.mode, abs.method = abs.method)
   if(absolute.mode) message(paste0("\n", ">>> Running ", "SVM in absolute mode"))

   eset<-as.data.frame(eset)
   # the authors reccomend to disable quantile normalizeation for RNA seq.
   # (see CIBERSORT website).
   quantile_norm = arrays
   ##############################
   res<-CIBERSORT(sig_matrix = reference,
                  mixture_file = eset,
                  perm = perm,
                  QN = quantile_norm,
                  absolute = absolute.mode,
                  abs_method = abs.method)
   ###############################

  }else if(method=="lsei"){

    message(paste0("\n", ">>> Running ", "cell estimation in lsei mode"))
    eset <- as.matrix(eset)
    reference <- as.matrix(reference)
    # quantile normalized
    if (arrays){
      roweset <- rownames(eset)
      coleset <- colnames(eset)
      eset <- normalize.quantiles(eset)
      rownames(eset) <- roweset
      colnames(eset) <- coleset
    }
    # scale_reference
    if (scale_reference){
      reference <- (reference - mean(reference))/sd(as.vector(reference))
    }
    Ymedian <- max(median(eset),1)
    # common eset
    common <- intersect(rownames(eset), rownames(reference))
    eset <- eset[match(common, rownames(eset)), ]
    reference <- reference[match(common, rownames(reference)), ]
    # deconvolution
    output <- matrix()
    # message(paste0("\n", ">>> Running ", "cell estimation in lsei mode"))
    Numofx <- ncol(reference)
    AA <- reference
    EE <- rep(1, Numofx)
    FF <- 1
    GG <- diag(nrow=Numofx)
    HH <- rep(0, Numofx)
    out.all <- c()
    itor <- 1
    samples <- ncol(eset)
    while (itor <= samples){
        BB <- eset[, itor]
        BB <- (BB - mean(BB))/sd(BB)
        out <- lsei(AA, BB, EE, FF, GG, HH)
        out.all <- rbind(out.all, out$X)
        itor <- itor + 1
      }
    rownames(out.all) <- colnames(eset)
    res <- out.all
  }

  ###############################
  colnames(res)<-gsub(colnames(res),pattern = "\\.",replacement = "\\_")
  colnames(res)<-gsub(colnames(res),pattern = "\\ ",replacement = "\\_")
  colnames(res)<-paste0(colnames(res),"_",method)
  res<-as.data.frame(res)

  if(!is.null(project)){
    res$ProjectID<-project
    res<-res[,c(ncol(res),1:ncol(res)-1)]
  }

  res<-tibble::rownames_to_column(res,var = "ID")
  return(res)

}


#' Deconvoluting using the TIMER technique
#' Unlike the other methods, TIMER needs the specification of the cancer type for each sample.
#'
#' @param eset  gene matrix
#' @param project default is NULL, project name used to distinguish different data sets
#' @param indications a n-vector giving and indication string (e.g. 'brca') for each sample. Accepted indications are 'kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg','lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct','ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca','uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol'
#'
#' @return
#' @export
#'
#' @examples
#'
deconvo_timer = function(eset,project = NULL,indications = NULL) {

  indications = tolower(indications)
  checkmate:: assert("indications fit to mixture matrix", length(indications) == ncol(eset))
  args = new.env()
  args$outdir = tempdir()
  args$batch = tempfile()
  lapply(unique(indications), function(ind) {
    tmp_file = tempfile()
    tmp_mat = eset[, indications == ind, drop=FALSE] %>% as_tibble(rownames = "gene_symbol")
    readr:: write_tsv(tmp_mat, tmp_file)
    cat(paste0(tmp_file, ",", ind, "\n"), file=args$batch, append=TRUE)
  })
  # reorder results to be consistent with input matrix
  results <- deconvolute_timer.default(args)[, make.names(colnames(eset))]


  colnames(results) <- colnames(eset)
  results<-as.data.frame(t(results))
  colnames(results)<-paste(colnames(results),"_TIMER",sep = "")
  colnames(results)<-gsub(colnames(results),pattern = "\\.",replacement = "\\_")
  colnames(results)<-gsub(colnames(results),pattern = "\\ ",replacement = "\\_")

  if(!is.null(project)){
    results$project<-project
    results<-results[,c(ncol(results),1:ncol(results)-1)]
  }

  results<-rownames_to_column(results,var = "ID")

  return(results)
}



#' Deconvoluting micrornvironment using the quanTIseq technique
#'
#' @param eset  gene expression data
#' @param tumor logistic variable
#' @param arrays logistic variable, is data generated from microarray
#' @param scale_mrna logistic variable
#' @param project default is NULL, project name used to distinguish different data sets
#'
#' @return
#' @export
#'
#' @examples
deconvo_quantiseq = function(eset, project = NULL, tumor, arrays, scale_mrna) {


  res = deconvolute_quantiseq.default(mix.mat = eset, tumor=tumor, arrays=arrays, mRNAscale = scale_mrna)

  res<-as.data.frame(res)
  rownames(res)<-NULL
  res<-column_to_rownames(res,var = "Sample")

  colnames(res)<-paste(colnames(res),"_quantiseq",sep = "")

  colnames(res)<-gsub(colnames(res),pattern = "\\.",replacement = "\\_")
  colnames(res)<-gsub(colnames(res),pattern = "\\ ",replacement = "\\_")

  if(!is.null(project)){
    res$ProjectID<-project
    res<-res[,c(ncol(res),1:ncol(res)-1)]
  }



  res<-rownames_to_column(res,var = "ID")

  return(res)
}



#' Deconvoluting Tumor microenvironment on a transcriptomic dataset
#'
#' @param eset A gene expression matrix
#'   Either: A numeric matrix or data.frame with HGNC gene symbols as rownames and sample identifiers as colnames. In both cases, data must be on non-log scale.
#' @param project project name used to distinguish different data sets, default is NULL
#' @param method a string specifying the method.
#' Supported methods are `mcp_counter`, `epic`, `xcell`, `cibersort`, `cibersort_abs`, `ips`, `quantiseq`, `estimate`,`timer`, `svm`,`lsei`，`timer`, `quantiseq`.
#' @param tumor logical. use a signature matrix/procedure optimized for tumor samples,
#'   if supported by the method. Currently affects `EPIC`
#' @param arrays Runs methods in a mode optimized for microarray data.
#'   Currently affects `CIBERSORT`, `svm` and `xCell`.
#' @param perm  set permutations for statistical analysis (≥100 permutations recommended).
#' Currently affects `CIBERSORT` and `svm_ref`
#' @param reference immune cell gene matrix; eg lm22, lm6 or can be generate using generateRef/generateRef_rnaseq
#' @param scale_reference a logical value indicating whether the reference be scaled or not. If TRUE, the value in reference file will be centered and scaled in row direction. Currently affects `svm` and `lsei` method
#' @param platform character string indicating platform type. Defaults to "affymetrix"
#' Currently affects `ESTIMATE` method
#' @param plot Currently affects `IPS` method
#' @param scale_mrna  logical. If FALSE, disable correction for mRNA content of different cell types.
#'   This is supported by methods that compute an absolute score (EPIC and quanTIseq)
#' @param group_list tumor type list of samples
#' @param absolute.mode Run CIBERSORT or SVM in absolute mode (default = FALSE)
#' @param abs.method if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#' @param ... arguments passed to the respective method
#'
#' @return `data.frame` with `ID` as first column and other column with the
#'     calculated cell fractions for each sample.
#' @author Dongqiang Zeng
#' @author Rongfang Shen
#' @name deconvo_tme
#' @export deconvo_tme
#' @examples
#'
deconvo_tme = function(eset,
                       project = NULL,
                       method = tme_deconvolution_methods,
                       arrays = FALSE,
                       tumor = TRUE,
                       perm = 1000,
                       reference,
                       scale_reference,
                       plot = FALSE,
                       scale_mrna,
                       group_list = NULL,
                       platform = "affymetrix",
                       absolute.mode = FALSE,
                       abs.method = "sig.score",
                       ...) {

  # message(paste0("\n", ">>> Running ", method))

  # run selected method
  res = switch(method,
               xcell = deconvo_xcell(eset, project ,arrays = arrays, ...),

               mcpcounter = deconvo_mcpcounter(eset, project, ...),

               epic = deconvo_epic(eset, project ,tumor = tumor, ...),

               cibersort = deconvo_cibersort(eset, project, absolute = FALSE, arrays = arrays,perm = perm, ...),

               cibersort_abs = deconvo_cibersort(eset, project, absolute = TRUE, abs_method = abs.method,arrays = arrays,perm = perm, ...),

               ips = deconvo_ips(eset,project,plot = plot, ...),

               quantiseq = deconvo_quantiseq(eset,project,tumor=tumor, arrays=arrays, scale_mrna=scale_mrna, ...),

               estimate = deconvo_estimate(eset,project,platform, ...),

               timer = deconvo_timer(eset,project, indications = group_list, ...),

               svm = deconvo_ref(eset, project, reference = reference, arrays = arrays, method = "svm", absolute.mode = absolute.mode,abs.method = abs.method,perm,...),

               lsei = deconvo_ref(eset, project, reference = reference, arrays = arrays, method = "lsei",scale_reference,perm,...) )

  res<-tibble::as_tibble(res)
  return(res)
}






