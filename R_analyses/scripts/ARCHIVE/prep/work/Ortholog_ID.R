library(GenomicFeatures)
library("readxl")
library(stringr)
library(dplyr)
library(tidyr)
library("AnnotationDbi")
library("org.Dm.eg.db")

dros_txdb_coords <- makeTxDbFromGFF("indata/ref_files/dros.gtf", format = 'gff')
k <- keys(dros_txdb_coords, keytype = "GENEID")
dros_txdf <- AnnotationDbi::select(dros_txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID") %>% 
  filter(is.na(TXNAME) == FALSE)
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
orths_t1 <- merge(orths_t0, TDel_txdf, by.x = 'TDel_TXN', by.y = 'TDel_TXN', all.y = T)
orths_t1 <- orths_t1[-which(duplicated(orths_t1$TDel_GID) & is.na(orths_t1$Dros_TXN)),]
products <- rtracklayer::import("indata/ref_files/TDel_GCF_002237135.1.gtf") %>% 
  as.data.frame() %>% 
  filter(gbkey == "CDS") %>% 
  dplyr::select(gene_id, product) %>% 
  distinct() %>% 
  filter(is.na(product) == FALSE) %>% 
  filter(duplicated(gene_id) == FALSE)

orths_t1$TDel_GID <- gsub("gene-", "", orths_t1$TDel_GID)
orths_t2 <- orths_t1 %>% 
  merge(products, by.x = 'TDel_GID', by.y = 'gene_id', all.x = TRUE, all.y = F)
orths_t2$product[is.na(orths_t2$product) == TRUE] <- "missing"


orths_t2 <- orths_t2 %>%
  mutate(product = str_replace_all(product, " ", "_"))

ortholog_table <- orths_t2 %>% 
  left_join(dros_txdf, by = 'Dros_TXN')

a <- ortholog_table %>% group_by(TDel_GID) %>% filter(n() > 1) %>% filter(is.na(Dros_GID) == F & 
                                                                        is.na(product) == F)
b <- ortholog_table %>% group_by(TDel_GID) %>% filter(n() == 1)
ortholog_table <- rbind(a,b)


save(ortholog_table, file = "outdata/RData/ortholog_table.RData")


##### NEXT THING TO DO IS LINK IN THE MARKER GENES 
markers <- read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx")
markers2 <- read_excel("indata/markers/flyatlas_dros_markers.xlsx") %>% 
  filter(Tissue == "testis")
markers3 <- readxl::read_excel("indata/markers/WITTPLOS2021/Compiled_markers13012022.xlsx", col_names = TRUE) %>% 
  pivot_longer(c("Marker1", "Marker2"), names_to = 'gene')
markers3 <- markers3[,c("Cell_type", "value")]
markers3 <- markers3[is.na(markers3$value) == FALSE,]

collapse <- function(x){
    gene <- t(str_split(x[15], ",", simplify = TRUE))
    gene <- gsub(" ", "", gene, fixed = TRUE)
    gene <- c(str_split(gene, ":", simplify = TRUE))
    gene <- gene[gene != ""]
    Cluster = rep(x[2], length(gene))
    return(data.frame())
}

markers2list <- apply(markers2, 1, function(x)(return(
  data.frame(
    Cluster = rep(x[2]), 
    gene = t(str_split(x[15], ",", simplify = TRUE)) %>% 
      gsub(pattern = " ", replacement =  "")) 
  )))

markers2df <- bind_rows(markers2list)
markers2df$gene <- gsub(" ", "", markers2df$gene, fixed = TRUE)
colnames(markers) <- paste("testis-", colnames(markers), sep = "")

ortholog_table2 <- merge(ortholog_table, markers, by.x = 'Dros_GID', by.y = 'testis-Gene', all.x = TRUE, 
                        )
colnames(markers2df) <- paste("atlas-", colnames(markers2df), sep = "")
ortholog_table3 <- merge(ortholog_table2, markers2df, by.x = 'Dros_GID', by.y = 'atlas-gene', all.x = TRUE)


save(ortholog_table1, ortholog_table2, ortholog_table3, file = "outdata/RData/orthologs.RData")


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

ortholog_table4 <- ortholog_table3 %>% 
  left_join(cell_cycle, by = c("Dros_GID" = "symbol"), keep = TRUE) 

cell_cycle_fin <- ortholog_table %>% 
  left_join(cell_cycle, by = c("Dros_GID" = "symbol")) %>% 
  filter(., is.na(phase) == FALSE)
write.table(cell_cycle_fin, file = "indata/markers/cell_cycle_markers_complete.csv")

atlas_clusters <- unique(ortholog_table4$`atlas-Cluster`)

ortholog_table4$comp_clusters <- ortholog_table4$`testis-Cluster`
ortholog_table4$comp_clusters[which(ortholog_table4$`atlas-Cluster` %in% 
                                      atlas_clusters[grep('spermatoc', unique(ortholog_table4$`atlas-Cluster`))] == T)] <- "Spermatocytes"
ortholog_table4$comp_clusters[which(ortholog_table4$`atlas-Cluster` %in% 
                               atlas_clusters[grep('cyst', unique(ortholog_table4$`atlas-Cluster`))] == T)] <- "Cyst"
ortholog_table4$comp_clusters[which(ortholog_table4$`atlas-Cluster` %in% 
                                     atlas_clusters[grep('epith', unique(ortholog_table4$`atlas-Cluster`))] == T)] <- "Epithelial cells"
ortholog_table4$comp_clusters[which(ortholog_table4$`atlas-Cluster` %in% 
                                      atlas_clusters[grep('hub', unique(ortholog_table4$`atlas-Cluster`))] == T)] <- "Hub cells"
ortholog_table4$comp_clusters[which(ortholog_table4$`atlas-Cluster` %in% 
                                      atlas_clusters[grep('spermatid', unique(ortholog_table4$`atlas-Cluster`))] == T)] <- "Spermatids"
ortholog_table4$comp_clusters[which(ortholog_table4$`atlas-Cluster` %in% 
                                      atlas_clusters[grep('spermatogon', unique(ortholog_table4$`atlas-Cluster`))] == T)] <- "GSC, Early spermatogonia"


ortholog_table4$comp_clusters[which(ortholog_table4$`testis-Cluster` %in% 
                                      c("Early spermatids","Mature spermatids") == T)] <- "Spermatids"
ortholog_table4$comp_clusters[which(ortholog_table4$`testis-Cluster` %in% 
                                      c("GSC, Early spermatogonia","Late spermatogonia") == T)] <- "GSC, Early spermatogonia"
ortholog_table4$comp_clusters[which(ortholog_table4$`testis-Cluster` %in% 
                                      c("Early spermatocytes","Late spermatocytes") == T)] <- "Spermatocytes"








save(ortholog_table1, ortholog_table2, ortholog_table3, ortholog_table4, file = "outdata/RData/orthologs.RData")
write.table(ortholog_table4, "data/ortholog_table.txt")


