


#' Source code for the TIMER deconvolution method.
#'
#' This code is adapted from https://github.com/hanfeisun/TIMER, which
#' again is an adapted version of the original TIMER source code
#' from http://cistrome.org/TIMER/download.html.
#'
#' The method is described in Li et al. Genome Biology 2016;17(1):174., PMID 27549193.


#' TimerINFO
#'
#' @param string
#'
#' @return
#' @export
#'
#' @examples
TimerINFO <- function(string) {
  message(sprintf('## %s\n', string))
}




#' TIMER signatures are cancer specific. This is the list of available cancer types.
#'
#' @export
timer_available_cancers <- c('kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg',
                         'lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct',
                         'ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca',
                         'uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol')


#' Remove batch effect of expression set
#'
#' @param cancer.exp
#' @param immune.exp
#' @param immune.cellType
#'
#' @return
#' @export
#'
#' @examples
RemoveBatchEffect <- function(cancer.exp, immune.exp, immune.cellType) {
  ## intersect the gene names of cancer.exp and immune.exp
  tmp.dd <- as.matrix(cancer.exp)
  tmp <- sapply(strsplit(rownames(cancer.exp), '\\|'),
                function(x) x[[1]])
  rownames(tmp.dd) <- tmp
  tmp.dd <- as.matrix(tmp.dd[which(nchar(tmp)>1), ])
  tmp.ss <- intersect(rownames(tmp.dd), rownames(immune.exp))

  ## bind cancer and immune expression data into one dataframe
  N1 <- ncol(tmp.dd)

  tmp.dd <- cbind(tmp.dd[tmp.ss, ], immune.exp[tmp.ss, ])
  tmp.dd <- as.matrix(tmp.dd)
  mode(tmp.dd) <- 'numeric'

  ## remove batch effects
  N2 <- ncol(immune.exp)
  tmp.batch <- c(rep(1, N1), rep(2, N2))
  tmp.dd0 <-sva:: ComBat(tmp.dd, tmp.batch, c(), BPPARAM = BiocParallel::bpparam("SerialParam"))

  ## separate cancer and immune expression data after batch effect removing
  dd.br <- tmp.dd0[, 1:N1]
  immune.exp.br <- tmp.dd0[, (N1+1):(N1+N2)]

  ## a immune category has multiple samples, use the median expression level for a gene
  tmp0 <- c()
  for(kk in unique(names(immune.cellType))){
    tmp.vv <- which(names(immune.cellType)==kk)
    tmp0 <- cbind(tmp0, apply(immune.exp.br[, tmp.vv], 1, median, na.rm=T))
  }


  immune.exp.agg.br <- tmp0
  colnames(immune.exp.agg.br) <- unique(names(immune.cellType))
  return(list(as.matrix(dd.br), immune.exp.br, immune.exp.agg.br))
}




#' process batch table and check cancer types.
#'
#' @param args environment
#'
#' @return
#' @export
#'
#' @examples
check_cancer_types <- function(args) {
  if (length(args$batch) != 0) {
    TimerINFO("Enter batch mode\n")
    cancers <- as.matrix(read.table(args$batch, sep=","))
  } else {
    cancers<- c(args$expression, args$category)
    dim(cancers) <- c(1, 2)
  }
  # print(cancers)
  for (i in seq(nrow(cancers))) {
    cancer.category <- cancers[i, 2]
    if (!(cancer.category %in% timer_available_cancers)) {
      stop(paste('unknown cancers:', cancer.category))
    }
  }
  return(cancers)
}




#' Constrained regression method implemented in Abbas et al., 2009
#'
#' @param XX immume expression data
#' @param YY cancer expression data
#' @param w
#'
#' @return
#' @export
#'
#' @examples
GetFractions.Abbas <- function(XX, YY, w=NA){
  ss.remove=c()
  ss.names=colnames(XX)
  while(T){
    if(length(ss.remove)==0)tmp.XX=XX else{
      if(is.null(ncol(tmp.XX)))return(rep(0, ncol(XX)))
      tmp.XX=tmp.XX[, -ss.remove]
    }
    if(length(ss.remove)>0){
      ss.names=ss.names[-ss.remove]
      if(length(ss.names)==0)return(rep(0, ncol(XX)))
    }
    if(is.na(w[1]))tmp=lsfit(tmp.XX, YY, intercept=F) else tmp=lsfit(tmp.XX, YY, w, intercept=F)
    if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
    if(min(tmp.beta>0))break
    ss.remove=which.min(tmp.beta)
  }
  tmp.F=rep(0, ncol(XX))
  names(tmp.F)=colnames(XX)
  tmp.F[ss.names]=tmp.beta
  return(tmp.F)
}

#' Convert Rowname To Loci
#'
#' @param cancerGeneExpression
#'
#' @return
#' @export
#'
#' @examples
ConvertRownameToLoci <- function(cancerGeneExpression) {
  ## Extract only the loci information for row name

  ## Example of origin row name is 'LOC389332|389332'
  ## Coverted row name is 'LOC389332'

  ## Args:
  ##   geneExpression: the orginal geneExpression load from .Rdata file
  ##
  ## Returns:
  ##   Modified geneExpression

  tmp <- strsplit(rownames(cancerGeneExpression), '\\|')
  tmp <- sapply(tmp, function(x) x[[1]])
  tmp.vv <- which(nchar(tmp) > 1)
  rownames(cancerGeneExpression) <- tmp
  extracted <- as.matrix(cancerGeneExpression[tmp.vv, ])
  colnames(extracted) <- colnames(cancerGeneExpression)
  return(extracted)
}



#' Input gene expression
#'
#' @param path path of data
#'
#' @return
#' @export
#'
#' @examples
ParseInputExpression <- function(path) {
  ret <- read.csv(path, sep='\t', row.names=1)
  ret <- as.matrix(ret)
  mode(ret) <- 'numeric'
  # ret <- ConvertRownameToLoci(ret)
  return(ret)
}


#' Draw QQ Plot
#'
#' @param cancer.exp
#' @param immune.exp
#' @param name
#'
#' @return
#' @export
#'
#' @examples
DrawQQPlot <- function(cancer.exp, immune.exp, name='') {
  ## q-q plot by sample should look like a straight line.
  ## Extreme values may saturate for Affy array data, but most of the data should align well.
  qq <- qqplot(cancer.exp, immune.exp, xlab='Tumor Expression', ylab='Ref Expression',
         main='Sample-Sample Q-Q plot')
  mtext(name, col="gray11")

  # get part of the points for fit linear, remove bottom 40%, and top 10%
  start <- 0.4 * length(qq$x)
  end <- 0.9 * length(qq$x)
  qq.sub <- list(x=qq$x[start:end], y=qq$y[start:end])
  fit <-lm(y ~ x, data=qq.sub)
  abline(fit, col="blue")
}

#' Get Outlier Genes
#'
#' @param cancers
#'
#' @return
#' @export
#'
#' @examples
GetOutlierGenes <- function (cancers) {
  ## Return a union of  outlier genes.
  ## The top 5 expressed genes in each sample is treated as outlier here.
  outlier.total <- c()
  for (i in seq(nrow(cancers))) {
    cancer.expFile <- cancers[i, 1]
    cancer.expression <- ParseInputExpression(cancer.expFile)
    for (j in 1:ncol(cancer.expression)) {
      outlier <- rownames(cancer.expression)[tail(order(cancer.expression[,j]), 5)]
      outlier.total <- c(outlier.total, outlier)
    }
  }
  return(unique(outlier.total))
}



#' deconvolute tumor microenvironment using timer
#'
#' @param args environment
#'
#' @return
#' @export
#'
#' @examples
deconvolute_timer.default = function(args) {


  # data("cancer_type_genes")
  # data("immuneCuratedData")

  cancers = check_cancer_types(args)

  TimerINFO('Loading immune gene expression')

  immune <- immuneCuratedData
  immune.geneExpression <- immune$genes
  immune.cellTypes <- immune$celltypes

  # message(immune$celltypes)

  outlier.genes <- sort(GetOutlierGenes(cancers))


  print(paste("Outlier genes:", paste(outlier.genes, collapse=' ')))

  dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(paste(args$outdir, '/results', sep=''))) {
    dir.create(paste(args$outdir, '/results', sep=''))
  }

  abundance.score.matrix <- c()
  pdf(paste(args$outdir, '/results/output.pdf', sep=''))
  for (i in 1:nrow(cancers)) {
    cancer.expFile <- cancers[i, 1]
    cancer.category <- cancers[i, 2]
    # gene.selected.marker.path <- system.file("extdata", "timer", "precalculated", paste0("genes_", cancer.category, ".RData"),
    #                                          package = "IOBR", mustWork = TRUE)
    cancer.expression <- ParseInputExpression(cancer.expFile)
    index <- !(row.names(cancer.expression) %in% outlier.genes)
    cancer.expression <- cancer.expression[index, , drop=FALSE]
    cancer.colnames <- colnames(cancer.expression)

    TimerINFO(paste("Removing the batch effect of", cancer.expFile))
    for (j in 1:length(cancer.colnames)) {
      DrawQQPlot(cancer.expression[,j], immune.geneExpression[,1], name=cancer.colnames[j])
    }

    tmp <- RemoveBatchEffect(cancer.expression, immune.geneExpression, immune.cellTypes)
    cancer.expNorm <- tmp[[1]]
    immune.expNormMedian <- tmp[[3]]

    for (j in 1:length(cancer.colnames)) {
      DrawQQPlot(cancer.expNorm[,j], immune.expNormMedian[,1],
               name=paste("After batch removing and aggregating for", cancer.colnames[j]))
    }

    gene.selected.marker <- cancer_type_genes[[which(names(cancer_type_genes)==cancer.category)]]
    gene.selected.marker <- intersect(gene.selected.marker, row.names(cancer.expNorm))
    XX = immune.expNormMedian[gene.selected.marker, c(-4)]
    YY = cancer.expNorm[gene.selected.marker, , drop=FALSE]

    for (j in 1:length(cancer.colnames)) {
      fractions <- GetFractions.Abbas(XX, YY[,j])
      # print (paste("Fractions for", cancer.expFile, cancer.colnames[j]))
      # print (fractions)
      barplot(fractions, cex.names=0.8, names.arg=names(fractions), xlab="cell type", ylab="abundance",
          main=paste("Abundance estimation for", cancer.colnames[j]))
      box()

      abundance.score.matrix <- cbind(abundance.score.matrix, fractions)
      colnames(abundance.score.matrix)[ncol(abundance.score.matrix)] <- cancer.colnames[j]
    }

  }

  dev.off()
  write.table(abundance.score.matrix, paste(args$outdir, '/results/score_matrix.txt', sep=''),
      sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

  return(abundance.score.matrix)
}

