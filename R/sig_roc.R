






#' Drawing multiple AUC into one graph
#'
#' @param data data frame with response and variables
#' @param variables variables
#' @param fig.path default is current working directory
#' @param palette default is `jama`, if number of variables is bigger than colors of palettes, `random` color will be applied
#' @param cols users can provide cols manually
#' @param main main title of plot
#' @param alpha default is 1
#' @param response response name of data
#' @param file.name file.name, default is
#' @param compare default is FALSE
#' @param compare_method  default is bootstrap, other option: “delong”, “venkatraman”
#' @param boot.n default is 100 when bootstrap is chosen.
#' @param smooth default is TRUE
#'
#' @author Dongqiang Zeng
#' @return
#' @export
#'
#' @examples
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
                       legacy.axes=T)

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
                  legacy.axes=T)

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
