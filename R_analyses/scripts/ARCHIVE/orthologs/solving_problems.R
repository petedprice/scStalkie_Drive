library(GenomicFeatures)
library("readxl")
library(stringr)
library(dplyr)
library(tidyr)
dros_txdb_coords <- makeTxDbFromGFF("indata/ref_files/dros.gtf", format = 'gff')
k <- keys(dros_txdb_coords, keytype = "GENEID")
dros_txdf <- AnnotationDbi::select(dros_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM", "TXSTART", "TXEND"), keytype = "GENEID")
colnames(dros_txdf) <- c("Dros_GID", "Dros_TXN", "Dros_CHR", "start", "end")
dros_txdf$length <- abs(dros_txdf$start - dros_txdf$end)
markers <- readxl::read_excel("indata/markers/WITTPLOS2021/Compiled_markers13012022.xlsx")
goi <- unique(unlist(c(markers[,2:ncol(markers)])))

ss <- filter(dros_txdf, Dros_GID %in% goi)
fileConn <- file("indata/markers/del.txt")
writeLines((c(ss$Dros_TXN)), fileConn)
close(fileConn)

write.table(c(ss$Dros_TXN), "indata/markers/del.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "")
length(unique(ss$Dros_GID))

tois <- read.table("indata/markers/problems/toi.txt")
tois2 <- c("NM_001274464.1", "NM_001103704.3","NM_001259765.2", 
          "NM_134310.3", 'NM_057386.5', "NM_001273529.2",
          "NM_057285.4", "NM_142101.3","NM_132988.3","NM_170060.2",
          "NM_001169247.1", "NM_057452.4")


TDel_txdb_coords <- makeTxDbFromGFF("indata/ref_files/stalkie.gtf", format = 'gtf')
k <- keys(TDel_txdb_coords, keytype = "GENEID")
TDel_txdf <- AnnotationDbi::select(TDel_txdb_coords, keys = k,  columns = 
                                     c("TXNAME", "TXCHROM", "TXSTART", "TXEND"), keytype = "GENEID")
colnames(TDel_txdf) <- c("TDel_GID", "TDel_TXN", "TDel_CHR", "start", "end")

filter(orthologs_dros_atlas, (tois2 %in% Dros_TXN))

check_these_ones <- tois2[!(tois2 %in% orthologs_dros_atlas$Dros_TXN)]
  

high_matching_hits <- c("rna-XM_038095029.1",  "rna-XM_038096840.1")
filter(TDel_txdf, TDel_TXN %in% high_matching_hits)
nrow(TDel_txdf)
