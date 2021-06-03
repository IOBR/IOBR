



#' remove batch effect of two or three expression set
#'
#' @param eset1 expression set-1
#' @param eset2 expression set-2
#' @param eset3 expression set-3, default is null
#' @param labels labels for expression sets
#' @param check_eset default is true
#' @param palette palette from IOBR::palettes, default is nature
#'
#' @return
#' @export
#'
#' @examples
remove_batcheffect<-function(eset1, eset2, eset3 = NULL, check_eset = TRUE, labels = NULL, palette = "nrc"){


  cols<-IOBR::palettes(category = "box", palette = palette, show_col = F, show_message = F)

  eset1<-log2eset(eset1)
  if(check_eset) check_eset(eset1)

  eset2<-log2eset(eset2)
  if(check_eset) check_eset(eset2)


  if(!is.null(eset3)) {
    eset3<-log2eset(eset3)
    if(check_eset) check_eset(eset3)
    cols<-cols[1:3]
  }else{
    cols<-cols[1:2]
  }

  if(is.null(eset3)){
    comgene <- intersect(rownames(eset1), rownames(eset2))
    combined.expr <- cbind.data.frame(eset1[comgene,],
                                      eset2[comgene,])
    batch <- data.frame(batch = rep(c("eset1","eset2"), times = c(ncol(eset1),ncol(eset2))))
  }
  ######################

  if(!is.null(eset3)){
    comgene <- intersect(intersect(rownames(eset1), rownames(eset2)), rownames(eset3))
    combined.expr <- cbind.data.frame(eset1[comgene,],
                                      eset2[comgene,],
                                      eset3[comgene,])
    batch <- data.frame(batch = rep(c("eset1","eset2","eset3"), times = c(ncol(eset1),ncol(eset2),ncol(eset3))))

  }

  ######################

  # outfile = file.path(fig.dir, paste("PCA before and after batch effect manipulation", ".pdf",sep=""))
  # pdf(file=outfile, width = 7, height = 3.5)

  par(mfcol = c(1,2))

  if(is.null(labels) & is.null(eset3)){
    labels <-c("eset1", "eset2")
  } else if(is.null(labels) & !is.null(eset3)){
    labels <-c("eset1", "eset2", "eset3")
  }

  if(!is.null(eset3)){
    batch_for_pca<-rep(labels, times = c(ncol(eset1),ncol(eset2),ncol(eset3)))
  }else{
    batch_for_pca<-rep(labels, times = c(ncol(eset1),ncol(eset2)))
  }

  batchPCA(indata = t(scale(t(combined.expr))),
           batch = batch_for_pca,
           main.title = "Raw PCA for combined expression profile",
           cols =  cols,
           showID = F,
           cex = 0.7,
           showLegend = T)


  modcombat = model.matrix(~1, data = batch)
  combined.expr.combat <- as.data.frame(sva:: ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))

  batchPCA(indata = t(scale(t(combined.expr.combat))),
           batch = batch_for_pca,
           main.title = "Combat PCA for combined expression profile",
           cols =  cols,
           showID = F,
           cex = 0.7,
           showLegend = T)

  # invisible(dev.off())

  return(combined.expr.combat)
}
