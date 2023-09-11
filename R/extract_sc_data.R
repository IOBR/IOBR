






#' Extract data frame
#'
#' Extract and combine data frame with cells as rows and features as columns from Seurat assay data.
#' @param sce Seurat object
#' @param vars Vector of variables as features
#' @param assay Assay to pull data from
#' @param slot Specific assay data to get or set
#' @param combine_meta_data Whether to combine metadata with the extracted data frame
#'
#' @return Data frame with cells as rows and features as columns
#' @export
#'
#' @author Dongqiang Zeng
#' @examples
#' #Load data
#' data("pbmc_small")
#' pbmc_small
#' #Choose features
#' vars<-c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3", "PF4", "MYL9", "GNLY", "TREML1", "CA2")
#' #Get extracted data frame
#' eset<-extract_sc_data(sce= pbmc_small, vars= vars, assay= "RNA")
extract_sc_data<-function(sce, vars = NULL, assay, slot = "scale.data", combine_meta_data = TRUE){


  exist<-Seurat::Assays(sce)
  message(paste0(">>>--- Assays of seurat object: "))
  message(paste(">>>---", exist, collapse = " "))
  assay<-assay[assay%in%exist]

  if(length(assay)==0) stop(">>>--- There no assay in object...")

  eset_cbind<-data.frame("ID" = rownames(sce@meta.data), "index" = 1:nrow(sce@meta.data))

  for(i in 1:length(assay)){

    method<-assay[i]
    DefaultAssay(sce)<-method

    eset<- SeuratObject:: GetAssayData(sce, assay = assay, slot = slot)


    feas <- rownames(eset)[rownames(eset)%in% unique(vars)]
    if(length(feas)==0) stop(">>>== The required variables are not present in the expression matrix")
    eset<- eset[feas, ]
    # print(head(eset))
    eset<- as.data.frame(t(as.matrix(eset)))

    if(length(vars)==1) {

      eset <-as.data.frame(eset)
      rownames(eset) <- vars
      eset <- t(eset)
      eset <-as.data.frame(eset)
    }

    # print(head(eset[, 1:5]))
    eset<- tibble:: rownames_to_column(eset, var = "ID")

    # base::as.data.frame() %>%
    # tibble:: rownames_to_column(.,var = "id") %>%
    # dplyr:: filter(.$id%in%unique(vars)) %>%
    # column_to_rownames(.,var = "id") %>%
    # base:: t() %>%
    # base:: as.data.frame() %>%

    exit_vars<-vars[vars%in%colnames(eset)]

    if(length(exit_vars)==0) next

    #print(head(eset))

    if(length(assay)>1) colnames(eset)[2:ncol(eset)] <-paste0(colnames(eset)[2:ncol(eset)], "_", method)

    # colnames(eset)<-gsub(colnames(eset), pattern = "-", replacement =  "_")
    # colnames(eset)<-gsub(colnames(eset), pattern = "/", replacement = "_")
    # colnames(eset)<-gsub(colnames(eset), pattern = " ", replacement = "_")
    eset<-as.data.frame(eset)

    if(length(assay)==1){
      eset_cbind<-eset
    }else{
      eset_cbind<-inner_join(eset_cbind, eset, by ="ID")
    }

  }

  if(combine_meta_data){

    message(paste0(">>>--- Merging metadata... "))
    meta.data<-rownames_to_column(sce@meta.data, var = "ID")
    eset_cbind<- inner_join(meta.data, eset_cbind, by = "ID" )
  }

  # eset_cbind<-eset_cbind[,-which(colnames(eset_cbind)=="index")]
  #print(head(eset_cbind))
  return(eset_cbind)

}
