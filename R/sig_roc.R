






#' Drawing multiple AUC into one graph
#' 
#' This function `sig_roc` creates ROC curves for evaluating the performance of one or more predictors in a binary classification task. It is primarily designed to assess how well variables can distinguish between two classes, such as disease vs. no disease. The function generates a plot for each variable and optionally compares their Areas Under the Curve (AUCs) using statistical tests.
#'
#' @param data A data frame containing the response variable and predictors.
#' @param response The name of the response variable in the data frame, typically a binary outcome.
#' @param variables A vector of column names from `data` representing predictors for which ROC curves will be plotted.
#' @param fig.path Optional; the path where the output plot will be saved, defaults to the current working directory.
#' @param palette default is `jama`, if number of variables is bigger than colors of palettes, `random` color will be applied
#' @param cols Optional; a vector of colors to use for the curves. If not provided, colors are taken from the specified palette.
#' @param main The main title of the plot.
#' @param alpha Transparency level of the plot colors; default is 1 (opaque).
#' @param file.name Optional; the base name for the output file, defaults to "0-ROC of multiple variables" if not specified.
#' @param compare Logical; if TRUE, compares the ROC curves using statistical tests; default is FALSE.
#' @param compare_method  The method used for comparison if `compare` is TRUE; defaults to "bootstrap". Other options include "delong" and "venkatraman".
#' @param boot.n Number of bootstrap replicates if the bootstrap method is used; default is 100.
#' @param smooth Logical; if TRUE, the ROC curves are smoothed; default is TRUE
#'
#' @author Dongqiang Zeng
#' @return A list containing the AUC values, legend names for the plot, and optionally the comparison results.
#' @export 
#'
#' @examples
#' data("tcga_stad_pdata", package = "IOBR")
#' sig_roc(data = tcga_stad_pdata, response = "OS_status", variables = c("TMEscore_plus", "GZMB", "GNLY"))
sig_roc<-function(data,
                  response,
                  variables,
                  fig.path = ".",
                  main = NULL,
                  file.name = NULL,
                  palette = "jama",
                  cols = NULL,
                  alpha = 1,
                  compare = FALSE,
                  smooth  = TRUE,
                  compare_method = "bootstrap",
                  boot.n = 100){
  options(pROCProgress = list(name = "none"))
  data<-as.data.frame(data)

  variables<-variables[variables%in%colnames(data)]

  input<-as.data.frame(data[,c(response,variables)])

  var_counts<-length(variables[variables%in%colnames(data)])

  if(is.null(cols)){
    cols<- palettes(palette = palette, alpha = alpha, show_col = FALSE, show_message = FALSE)

    if(var_counts> length(cols)){
      cols<- palettes(category = "random", alpha = alpha, show_col = FALSE, show_message = FALSE)
    }
  }

  ########################################

  auc.out <- c()

  if(is.null(file.name)) file.name<- "0-ROC of multiple variables"

  outfile = file.path(fig.path, paste(file.name, ".pdf",sep=""))

  pdf(file = outfile, width = 5, height = 5)
  x <- pROC:: plot.roc(input[,1],input[,2],ylim=c(0,1),xlim=c(1,0),
                       smooth= smooth,
                       ci=TRUE,
                       main = main ,
                       col=cols[2],
                       lwd=1.5,
                       legacy.axes=T,
                       xlab = "False Positive Rate",
                       ylab = "True Positive Rate")

  ci.lower <- round(as.numeric(x$ci[1]),3)
  ci.upper <- round(as.numeric(x$ci[3]),3)

  auc.ci <- c(colnames(input)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
  ##############################
  for (i in 3:ncol(input)){
    x <-pROC::plot.roc(input[,1],input[,i],
                  add        =T,
                  smooth     =smooth,
                  ci         =TRUE,
                  col        =cols[i],
                  lwd        =2,
                  legacy.axes=T,
                  xlab = "False Positive Rate",
                  ylab = "True Positive Rate")

    ci.lower <- round(as.numeric(x$ci[1]),3)
    ci.upper <- round(as.numeric(x$ci[3]),3)
    auc.ci <- c(colnames(input)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out <- rbind(auc.out,auc.ci)
  }

  auc.out <- as.data.frame(auc.out)
  colnames(auc.out) <- c("Name","AUC","AUC CI")

  ################################
  if(compare){

    p.out <- c()

    for (i in 2:(ncol(input)-1)){
      for (j in (i+1):ncol(input)){
        p <-pROC:: roc.test(input[,1],input[,i],input[,j], method = compare_method,boot.n = boot.n, progress = "none")
        p.tmp <- c(colnames(input)[i],colnames(input)[j],p$p.value)
        p.out <- rbind(p.out,p.tmp)
      }
    }
    p.out <- as.data.frame(p.out)
    # head(p.out)
    colnames(p.out) <- c("ROC1","ROC2","p.value")
  }

  legend.name <- paste(colnames(input)[2:length(input)]," AUC = ",auc.out$AUC,sep=" ")
  legend("bottomright",
         legend=legend.name,
         col = cols[2:length(input)],
         lwd = 2,
         bty="n")

  dev.off()

  if(compare){
    return(list(auc.out = auc.out, legend.name = legend.name, p.out = p.out ))
  }else{
    return(list(auc.out = auc.out, legend.name = legend.name ))
  }


}
