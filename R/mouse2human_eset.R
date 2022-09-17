




#' Converting muouse gene symbol to human gene symbol of expression set.
#'
#' @param eset expression matrix or data frame
#' @param source default is ensembl, if an error is reported, set parameter to `local` and mus_human_gene_symbol will be use to convert gene symbols
#' @param is_matrix default is T, if not, column_of_symbol must be se
#' @param column_of_symbol default is null
#'
#' @return
#' @export
#'
#' @examples
mouse2human_eset<-function(eset, source = "ensembl", is_matrix = T, column_of_symbol = NULL, verbose = FALSE){


  if(is_matrix){
    genes<-rownames(eset)
  }else{
    eset<-remove_duplicate_genes(eset = eset, column_of_symbol = column_of_symbol)
    genes<-rownames(eset)
  }


  if(source=="ensembl"){

    require("biomaRt")
    ensembl = biomaRt:: useEnsembl(biomart="ensembl")
    if(verbose)  print(head(listDatasets(ensembl)))
    # Basic function to convert mouse to human gene names

    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    probe_data = getLDS(attributes = c("mgi_symbol"),
                        filters = "mgi_symbol",
                        values = genes ,
                        mart = mouse,
                        attributesL = c("hgnc_symbol"),
                        martL = human,
                        uniqueRows=T)
    colnames(probe_data)<-c("gene_symbol_mus", "gene_symbol_human")

  }else if(source=="local"){

    probe_data<-mus_human_gene_symbol

  }

  eset<-anno_eset(eset = eset, annotation = probe_data, symbol = "gene_symbol_human", probe = "gene_symbol_mus",method = "mean")
  return(eset)
}
