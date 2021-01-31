



#' Choose colors.
#'
#' @param category choose what type of plot you want to color,There are seven categories you can choose: `box`, `continue2`, `continue`, `random`, `heatmap`, `heatmap3`, `tidyheatmap`
#' @param palette Each category has some options to choose
#' @param alpha intensity of the color
#' @param counts number of colors, only affect continue colors
#' @param show_col logical variable, if TRUE, color will be print and show in the R studio
#' @param show_message print palettes that can be choose
#'
#' @author Dongqiang Zeng
#' @return
#' @export
#'
#' @examples
palettes<-function(category = "box", palette = "nrc",alpha = 1,counts = 50, show_col = TRUE, show_message = FALSE){

  message(paste0("There are seven categories you can choose: box, continue2, continue, random, heatmap, heatmap3, tidyheatmap "))


  if(category == "box"){
     if(show_message) message(paste0("There are ten palettes you can choose: nrc, jama, aaas, jco, paired1-4, accent, set2"))
    if(palette=="nrc"){
      mypal = ggsci:: pal_npg("nrc", alpha = alpha)(9)

    }else if(palette == "jama"){
      mypal = ggsci:: pal_jama(palette = c("default"),alpha = alpha)(7)

    }else if(palette == "aaas"){
      mypal = ggsci::pal_aaas(palette = c("default"),alpha = alpha)(9)

    }else if(palette == "jco"){
      mypal = ggsci::pal_jco(palette = c("default"),alpha = alpha)(9)

    }else if(palette == "paired1"){
      mypal = RColorBrewer:: brewer.pal(11,"Paired")
      mypal<-mypal[1:8]

    }else if(palette == "paired2"){
      mypal = RColorBrewer:: brewer.pal(11,"Paired")
      mypal<-mypal[3:10]

    }else if(palette == "paired3"){
      mypal = RColorBrewer::brewer.pal(11,"Paired")
      mypal<-mypal[5:11]

    }else if(palette == "paired4"){
      mypal = RColorBrewer::brewer.pal(11,"Paired")
      mypal<-mypal[7:11]

    }else if(palette == "accent"){
      mypal =RColorBrewer:: brewer.pal(8,"Accent")

    }else if(palette == "set2"){
      mypal =RColorBrewer:: brewer.pal(7,"Set2")


    }

    if(show_col){
      print(paste0("'", mypal,"'",collapse = ", "))
      scales:: show_col(mypal)
    }
  }

  if(category == "continue2"){
    if(show_message)  message(paste0("There are five palettes you can choose: nrc, jama, aaas, jco, rdbu"))
    if(palette=="nrc"){
      mypal = ggsci:: pal_npg("nrc", alpha = alpha)(9)
      mypal<-mypal[c(4,1)]

    }else if(palette == "jama"){
      mypal = ggsci:: pal_jama(palette = c("default"),alpha = alpha)(9)
      mypal<-mypal[c(1,4)]

    }else if(palette == "aaas"){
      mypal = ggsci::pal_aaas(palette = c("default"),alpha = alpha)(9)
      mypal<-mypal[c(1,6)]

    }else if(palette == "jco"){
      mypal = ggsci::pal_jco(palette = c("default"),alpha = alpha)(9)
      mypal<-mypal[c(1,2)]

    }else if(palette == "rdbu"){
      mypal = RColorBrewer::brewer.pal(11,"RdBu")
      mypal<-mypal[c(10,2)]

    }
    if(show_col){
      print(paste0("'", mypal,"'",collapse = ", "))
      scales:: show_col(mypal)
    }

  }

  if(category=="random"){

    mypal<-c(RColorBrewer:: brewer.pal(8, "Accent"),
             RColorBrewer:: brewer.pal(8,"Dark2"),
               c("#224444","#e68a00","#33adff","#a6a6a6","#439373",
                 "#67001F","#1B9E77","#FB9A99","#FDBF6F","#92C5DE",
                 "#33A02C","#D1E5F0","#92C5DE","#4393C3",
                 "#B2DF8A","#CAB2D6","#56B4E9","#BC3C29FF"),
             RColorBrewer:: brewer.pal(8,"Set1"),
             RColorBrewer:: brewer.pal(8,"Set2"))
    if(palette==1){
      mypal<- mypal[1:15]
    }else if(palette == 2){
      mypal<- mypal[8:16]
    }else if(palette == 3){
      mypal<- mypal[10:20]
    }
    if(show_col){
      print(paste0("'", mypal,"'",collapse = ", "))
      scales:: show_col(mypal)
    }

  }

  if(category == "continue"){
    if(show_message)  message(paste0("There are four palettes you can choose: rdbu, puor, blues, reds"))
    if(palette == "rdbu" ){
    mypal<-  RColorBrewer::brewer.pal(11,"RdBu")

    }else if(palette == "puor"){
      mypal<- RColorBrewer:: brewer.pal(11,"PuOr")

    }else if(palette == "blues"){
      mypal<- RColorBrewer:: brewer.pal(11,"Blues")

    }else if(palette == "reds"){
      mypal<- RColorBrewer:: brewer.pal(11,"Reds")

    }
    if(show_col){
      print(paste0("'", mypal,"'",collapse = ", "))
      scales:: show_col(mypal)
    }
  }

  if(category == "heatmap"){

    if(show_message) message(paste0("There are five palettes you can choose: pheatmap, virids, blues, reds, peach"))

    if(palette == "pheatmap" ){
      mypal<- rev(colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(counts)

    }else if(palette == "peach"){
      mypal<-  colorRampPalette(c("#3182bd", "white", "#dd1c77"))(counts)

    }else if(palette == "blues"){
      mypal<- rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Blues")))(counts)

    }else if(palette == "virids"){
      mypal<- viridis:: inferno(counts)

    }else if(palette == "reds"){
      mypal<- rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Reds")))(counts)
    }
    if(show_col){
      print(paste0("'", mypal,"'",collapse = ", "))
      scales:: show_col(mypal)
    }

  }

  if(category=="heatmap3"){
    if(show_message) message(paste0("There are six palettes you can choose: pheatmap, virids, blues, reds, peach, normal"))

    if(palette == "pheatmap" ){
      mypal<- c("#4575B4","#FEF9B6","#D73027")


    }else if(palette == "peach"){
      mypal<-  c("#3182bd", "white", "#dd1c77")


    }else if(palette == "blues"){
      mypal<-  c("#F7FBFF", "88BEDC", "#084594")


    }else if(palette == "virids"){
      mypal<- c("#000004FF", "#CD4347FF", "#FCFFA4FF")


    }else if(palette == "reds"){
      mypal<- c("#FFF5F0", "FB7555", "#99000D")


    }else if(palette == "normal"){
      mypal<- c("navy","white","firebrick")

    }
    if(show_col){
      print(paste0("'", mypal,"'",collapse = ", "))
      scales:: show_col(mypal)
    }

  }

  if(category == "tidyheatmap"){

    if(show_message)  message(paste0("There are six palettes you can choose: 1, 2, 3, 4, 5, 6"))
    if(palette==1){
      mypal<-circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), viridis::magma(5))
    }else if(palette==2){
      mypal<-circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), rev(RColorBrewer:: brewer.pal(n = 5, name = "RdYlBu")))
    }else if(palette==3){
      mypal<-circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), rev(RColorBrewer::brewer.pal(n = 5, name = "RdYlGn")))
    }else if(palette==4){
      mypal<-circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3), rev(RColorBrewer::brewer.pal(n = 5, name = "Spectral")))
    }else if(palette==5){
      mypal<-circlize::colorRamp2(c(-3, -1.5, 0, 1.5, 3),  rev(RColorBrewer:: brewer.pal(n = 5, name = "PiYG")))
    }else if(palette==6){
      mypal<-circlize::colorRamp2(c(-3,  0,  3),  c("navy","white", "firebrick"))
    }

  }

  return(mypal)
}

