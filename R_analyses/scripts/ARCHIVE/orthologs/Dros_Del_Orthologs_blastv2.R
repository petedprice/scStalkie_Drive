library(GenomicFeatures)
library("readxl")
library(stringr)
library(dplyr)
library(tidyr)
library(flybaseR)

dros_txdb_coords <- makeTxDbFromGFF("indata/ref_files/dros.gtf", format = 'gff')
k <- keys(dros_txdb_coords, keytype = "GENEID")
dros_txdf <- AnnotationDbi::select(dros_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
colnames(dros_txdf) <- c("Dros_GID", "Dros_TXN", "Dros_CHR")



TDel_txdb_coords <- makeTxDbFromGFF("indata/ref_files/stalkie.gtf", format = 'gtf')
k <- keys(TDel_txdb_coords, keytype = "GENEID")
TDel_txdf <- AnnotationDbi::select(TDel_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
colnames(TDel_txdf) <- c("TDel_GID", "TDel_TXN", "TDel_CHR")

orths_t0 <- read.table("indata/orthologs/Jan2023_top_hits/top_hits2.txt", 
                       sep = ',', header = FALSE)

colnames(orths_t0) <- c("TDel_TXN", "Dros_TXN")
orths_t0[,1] <- str_split(orths_t0[,1], "lstalk_", simplify = TRUE)[,2]
orths_t0[,2] <- str_split(orths_t0[,2], "ldross_rna-", simplify = TRUE)[,2]
orths_t0[,2] <- str_split(orths_t0[,2], "[)]", simplify = TRUE)[,1]
orths_t1 <- merge(orths_t0, TDel_txdf, by.x = 'TDel_TXN', by.y = 'TDel_TXN')
ortholog_table <- merge(orths_t1, dros_txdf, by.x = 'Dros_TXN', 
                        by.y = 'Dros_TXN')


##### NEXT THING TO DO IS LINK IN THE MARKER GENES 
markers <- readxl::read_excel("indata/markers/WITTPLOS2021/Compiled_markers13012022.xlsx", col_names = TRUE) %>% 
  pivot_longer(c("Marker1", "Marker2", "Marker3"), names_to = 'gene')
markers <- markers[,c("Cell_type", "value")]
markers <- markers[is.na(markers$value) == FALSE,]
omarkers <- merge(markers, ortholog_table, by.x = "value", by.y= "Dros_GID")
save()

cell_cycle <- read.csv("indata/markers/cell_cycle_markers.csv")
library("AnnotationDbi")
library("org.Dm.eg.db")
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


