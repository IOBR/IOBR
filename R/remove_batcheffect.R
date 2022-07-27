



#' remove batch effect of two or three expression set
#'
#' @param eset1 expression set-1, normalized data
#' @param eset2 expression set-2, normalized data
#' @param eset3 expression set-3, default is null
#' @param labels labels for expression sets
#' @param check_eset default is true
#' @param palette palette from IOBR::palettes, default is nature
#' @param log2 default is TRUE
#' @param path default is NULL
#' @param save_plot default is FALSE, if true, path could be set to folder name
#'
#' @return
#' @export
#'
#' @examples
remove_batcheffect<-function(eset1, eset2, eset3 = NULL, log2 = TRUE, check_eset = TRUE, labels = NULL, palette = "nrc", path = NULL, save_plot = FALSE){


  cols<-IOBR::palettes(category = "box", palette = palette, show_col = F, show_message = F)

  if(is.null(path)){
    abspath<-NULL
  }else{
    if(!file.exists(path)) dir.create(path)
    abspath<-paste(getwd(),"/",file_store,"/",sep ="" )
  }

  if(log2) eset1<-log2eset(eset1)
  if(check_eset) check_eset(eset1)

  if(log2) eset2<-log2eset(eset2)
  if(check_eset) check_eset(eset2)


  if(!is.null(eset3)) {
    if(log2)  eset3<-log2eset(eset3)
    if(check_eset) check_eset(eset3)
    cols<-cols[1:3]
  }else{
    cols<-cols[1:2]
  }

  if(is.null(eset3)){
    comgene <- intersect(rownames(eset1), rownames(eset2))
    comgene<-comgene[!comgene==""]
    comgene<-comgene[!is.na(comgene)]

    combined.expr <- cbind.data.frame(eset1[comgene,],
                                      eset2[comgene,])
    batch <- data.frame(batch = rep(c("eset1","eset2"), times = c(ncol(eset1),ncol(eset2))))
  }
  ######################

  if(!is.null(eset3)){
    comgene <- intersect(intersect(rownames(eset1), rownames(eset2)), rownames(eset3))
    comgene<-comgene[!comgene==""]
    comgene<-comgene[!is.na(comgene)]
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
    table(batch_for_pca)
  }else{
    batch_for_pca<-rep(labels, times = c(ncol(eset1),ncol(eset2)))
    table(batch_for_pca)
  }

  batchPCA(indata = t(scale(t(combined.expr))),
           batch = batch_for_pca,
           main.title = "Raw PCA for combined expression profile",
           cols =  cols,
           showID = F,
           cex = 0.7,
           showLegend = T)


  modcombat = model.matrix(~ 1, data = batch)
  combined.expr.combat <- as.data.frame(sva:: ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))

  batchPCA(indata = t(scale(t(combined.expr.combat))),
           batch = batch_for_pca,
           main.title = "Combat PCA for combined expression profile",
           cols =  cols,
           showID = F,
           cex = 0.7,
           showLegend = T)

  # if(save_plot){
  #   pdf(paste0(abspath, "PCA-before-and-after-batch-effect.pdf"), width = 6.5, height = 3)
  # }

  # invisible(dev.off())

  return(combined.expr.combat)
}
