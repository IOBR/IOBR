


setwd("E:/18-Github/Organization")

###################################

my_gtf<- "E:/18-Github/Organization/gencode.vM32.annotation.gtf"
library(rtracklayer)
gtf <- import(paste0(my_gtf))
#################################################
gtf_df<- as.data.frame(gtf)
head(gtf_df)

table(gtf_df$type)
table(gtf_df$gene_type)

summary(duplicated(gtf_df$gene_id))
anno_gc_vm32<- gtf_df[!duplicated(gtf_df$gene_id), ]
anno_gc_vm32<- anno_gc_vm32[!is.na(anno_gc_vm32$gene_id), ]
#################################################

feas<-c("gene_id", "gene_name", "mgi_id", "gene_type", "start", "end", "transcript_id", "ont")
anno_gc_vm32<-anno_gc_vm32[,feas]
head(anno_gc_vm32)
colnames(anno_gc_vm32)<-c("id", "symbol", "mgi_id", "gene_type", "start", "end", "transcript_id", "ont")
save(anno_gc_vm32, file = "anno_gc_vm32.RData")
#################################################









# tested with R v4.0.0 and tidyverse v1.3.1
#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
read_gtf <- function(file) {

  require(tidyverse)
  cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")

  # read in raw gtf as tsv and remove comment rows
  messy <- read_tsv(file, col_names = F, comment = "#") %>%
    `colnames<-`(cnames)

  # get the unique attribute types
  # this assumes there are no spaces in the attribute names
  att_names <- messy %>%
    select(attribute) %>%
    apply(., MARGIN = 1, FUN = str_split, pattern = '"; ') %>%
    unlist() %>% trimws() %>% trimws(whitespace = ";") %>%
    sub(" .*$", "", .) %>% unique()

  att_names <- att_names[att_names != ""]

  # for each attribute type, create column
  # apply over gtf to fill in rows where attribute type is found
  for (att in att_names) {

    colatt <- apply(messy, MARGIN = 1, function(x) {

      var <- str_extract(string = x[9],
                         pattern = sprintf('";/s+%1$s[^;]+|^%1$s[^;]+;[^"]+"', att)) %>%
        trimws(whitespace = '["; ]+', which = 'left') %>%
        str_extract('(?<=")[^"]+(?=")')

    })

    messy <- messy %>% add_column("{att}" := colatt)

  }

  # remove original attribute column
  messy %>% select(-c(attribute))

}

