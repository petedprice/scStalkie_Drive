library(GenomicFeatures)
library("readxl")
library(stringr)

dros_txdb_coords <- makeTxDbFromGFF("indata/ref_files/dros.gtf", format = 'gff')
k <- keys(dros_txdb_coords, keytype = "GENEID")
dros_txdf <- AnnotationDbi::select(dros_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")

TDel_txdb_coords <- makeTxDbFromGFF("indata/ref_files/stalkie.gtf", format = 'gff')
k <- keys(TDel_txdb_coords, keytype = "GENEID")
TDel_txdf <- AnnotationDbi::select(TDel_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")


orthologs <- read.table("indata/orthologs/Orthofinder_oldref_practice_211122/Orthogroups/Orthogroups.tsv", sep = '\t', header = TRUE)


oto <- function(x){
  if (length(str_split(x[2], " ", simplify = TRUE)) < 2 & 
      length(str_split(x[3], " ", simplify = TRUE)) < 2){
    return(x)
  }
}
x <- orthologs[2,]
oto(orthologs[237,])
one_to_one <- unlist(apply(orthologs, 1, oto))
onto <- matrix(data = one_to_one, ncol = 3, byrow = TRUE) %>% as.data.frame()
colnames(onto) <- c("orthogroup", "Dros", "TDel")
onto[,2] <- str_split(onto[,2], "rna-", simplify = TRUE)[,2]
onto[,3] <- str_split(onto[,3], "rna-", simplify = TRUE)[,2]

Dr_TD <- merge(onto, dros_txdf[dros_txdf$TXNAME %in% onto$Dros,1:2], by.x = 'Dros', 
               by.y = 'TXNAME')

##### NEXT THING TO DO IS LINK IN THE MARKER GENES 
markers <- read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx")
Dr_TD_markers <- merge(Dr_TD, markers, by.x = 'GENEID', by.y = 'Gene', all.x = TRUE)

