library(GenomicFeatures)
library("readxl")
library(stringr)
library(dplyr)
library(tidyr)
library(flybaseR)
library(ggplot2)
library(tidyverse)

TDel_txdb_coords <- makeTxDbFromGFF("indata/ref_files/TDel_GCF_002237135.1.gtf", format = 'gtf')
k <- keys(TDel_txdb_coords, keytype = "GENEID")
TDel_txdf <- AnnotationDbi::select(TDel_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
colnames(TDel_txdf) <- c("TDel_GID", "TDel_TXN", "TDel_CHR")

orths_t0 <- read.table("indata/orthologs/Jan2023_top_hits/top_hits2.txt", 
                       sep = ',', header = FALSE)

a <- rtracklayer::import("indata/ref_files/TDel_GCF_002237135.1.gtf") %>% 
  as.data.frame()




colnames(orths_t0) <- c("TDel_TXN", "Dros_TXN")
orths_t0[,1] <- str_split(orths_t0[,1], "lstalk_rna-", simplify = TRUE)[,2]
orths_t0[,2] <- str_split(orths_t0[,2], "ldross_rna-", simplify = TRUE)[,2]
orths_t0[,2] <- str_split(orths_t0[,2], "[)]", simplify = TRUE)[,1]
ortholog_table <- left_join(orths_t0, TDel_txdf, by = c('TDel_TXN' = 'TDel_TXN')) %>% 
  left_join(dros_txdf, by = c('Dros_TXN' = 'Dros_TXN'))

ortholog_table %>% group_by(TDel_CHR) %>% 
  summarise(n = n())
ortholog_table %>% group_by(Dros_CHR) %>% 
  summarise(n = n())

ortholog_table %>% filter(TDel_CHR == 'NC_051846.1') %>% 
  ggplot(aes(x = Dros_CHR)) + geom_bar()

ortholog_table %>% filter(Dros_CHR == 'NC_004354.4') %>% 
  ggplot(aes(x = TDel_CHR)) + geom_bar()

NC_004354.4
##### NEXT THING TO DO IS LINK IN THE MARKER GENES 
markers <- readxl::read_excel("indata/markers/WITTPLOS2021/Compiled_markers13012022.xlsx", col_names = TRUE) %>% 
  pivot_longer(c("Marker1", "Marker2", "Marker3"), names_to = 'gene')
markers <- markers[,c("Cell_type", "value")]
markers <- markers[is.na(markers$value) == FALSE,]
omarkers <- merge(markers, ortholog_table, by.x = "value", by.y= "Dros_GID")
save


##### CHECKING WHICH GENES ARE ON THE X
DEL_txdb_coords <- makeTxDbFromGFF("../ref_files/T_dalmanni_v2_annotation.gff", format = 'gff')
k <- keys(DEL_txdb_coords, keytype = "GENEID")
del_txdf <- AnnotationDbi::select(dros_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
colnames(dros_txdf) <- c("Dros_GID", "Dros_TXN", "Dros_CHR")
dim(dros_txdf)

library("AnnotationDbi")
library("org.Dm.eg.db")


del <- read.table("../ref_files/T_dalmanni_v2_annotation.gff", sep = '\t')

dros <- read.table("indata/ref_files/dros.gtf", sep = '\t') %>% 
  subset(V3 == 'gene') 

dros$FB <- gsub(".*FLYBASE:", "", dros$V9) %>% 
  gsub(",GeneID.*", "", .) %>% 
  gsub(",.*", "", .)

del$FB <- str_split(del$V9, ";", simplify = T)[,2] %>% 
  gsub("Name=", "", .) %>% 
  gsub(",.*", "", .) %>% 
  gsub("[.].*", "", .)

dros <- dros[duplicated(dros$FB) == FALSE,] %>% 
  select(V1, FB) %>% 
  rename(CHR = V1)

del <- del[-duplicated(del$FB) == FALSE,] %>% 
  select(V1, FB) %>% 
  rename(CHR = V1)

dros_del <- dros %>% left_join(del, by = c("FB" = "FB"), suffix = c("dros", "del")) %>% 
  filter(is.na(CHRdel) == FALSE) 
dros_del[dros_del == 	'NC_004354.4'] <- "X"

dros_del %>% ggplot(aes(x = CHRdel, fill = CHRdros)) + geom_bar(position = 'dodge')
dros_del %>% ggplot(aes(fill = CHRdel, x = CHRdros)) + geom_bar(position = 'dodge')

del_X$symbol <- mapIds(org.Dm.eg.db, keys=del_X[,1], 
                       column="SYMBOL",  keytype="ENSEMBL",  multiVals="first") 

dros_X$symbol <- mapIds(org.Dm.eg.db, keys=dros_X[,1], 
                        column="SYMBOL",  keytype="ENSEMBL",  multiVals="first") 

del_X %>% filter(is.na(symbol) == FALSE)
dros_X %>% filter(is.na(symbol) == FALSE)


X_genes %>% filter(is.na(symbol) == FALSE) %>% 
  left_join(dros_txdf, by = c('symbol' == "Dros_GID"))
intersect(X_genes$symbol, dros_txdf$Dros_GID)



cell_cycle <- read.csv("indata/markers/cell_cycle_markers.csv")
cell_cycle$symbol <- mapIds(org.Dm.eg.db, keys=cell_cycle$geneID, 
                            column="SYMBOL",  keytype="ENSEMBL",  multiVals="first")

cell_cycle$entrez <- mapIds(org.Dm.eg.db,  keys=cell_cycle$geneID, 
                            column="ENTREZID",  keytype="ENSEMBL",  multiVals="first")

cell_cycle$name =   mapIds(org.Dm.eg.db, keys=cell_cycle$geneID, 
                           column="GENENAME", keytype="ENSEMBL", multiVals="first")

cell_cycle <- apply(cell_cycle,2,as.character) %>% as.data.frame()
cell_cycle <- ortholog_table %>% 
  left_join(cell_cycle, by = c("Dros_GID" = "symbol")) %>% 
  filter(., is.na(phase) == FALSE)
write.table(cell_cycle, file = "indata/markers/cell_cycle_markers_complete.csv")


