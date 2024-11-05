library(GenomicFeatures)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(ggpubr)

#### READING IN DATA -----
load("data/RData/testis.Cluster_marker_seurat.RData")

txdb_coords <- makeTxDbFromGFF("indata/ref_files/stalkie.gtf")
k <- keys(txdb_coords, keytype = "GENEID")
txdf <- AnnotationDbi::select(txdb_coords, keys = k,  columns = c("TXNAME", "TXCHROM"), keytype = "GENEID")
coords_genes <- as.data.frame(genes(txdb_coords, c("TXCHROM", "GENEID")))
coords_genes$GENEID <- str_split(coords_genes$GENEID, "gene-", simplify = T)[,2]
coords_genes$TXCHROM_name <- "unknown"
coords_genes$TXCHROM_name[coords_genes$TXCHROM == "NC_051846.1"] <- 1
coords_genes$TXCHROM_name[coords_genes$TXCHROM == "NC_051847.1"] <- 2
coords_genes$TXCHROM_name[coords_genes$TXCHROM == "NC_051848.1"] <- "X"
coords_genes$mid <- coords_genes$start + (coords_genes$width/2)




############# CHECKING distribution of expression on each chromosome across cell types and treatments 


Xgenes <- filter(coords_genes, TXCHROM_name == "X")$GENEID %>% intersect(rownames(seurat_marker))
genes1 <- filter(coords_genes, TXCHROM_name == 1)$GENEID %>% intersect(rownames(seurat_marker))
genes2 <- filter(coords_genes, TXCHROM_name == 2)$GENEID %>% intersect(rownames(seurat_marker))
genesUK <- filter(coords_genes, TXCHROM_name == 'unknown')$GENEID %>% intersect(rownames(seurat_marker))

total_genes <- nrow(seurat_marker)

DefaultAssay(seurat_marker) <- "RNA"
seurat_marker[['Xexp']] <- PercentageFeatureSet(object = seurat_marker, features = Xgenes)/
  (length(Xgenes) *100/total_genes)
seurat_marker[['CHR1exp']] <- PercentageFeatureSet(object = seurat_marker, features = genes1)/
  (length(genes1) *100/total_genes)
seurat_marker[['CHR2exp']] <- PercentageFeatureSet(object = seurat_marker, features = genes2)/
  (length(genes2) *100/total_genes)
seurat_marker[['UKexp']] <- PercentageFeatureSet(object = seurat_marker, features = genesUK)/
  (length(genesUK) *100/total_genes)

mean(seurat_marker$UKexp)
filter(seurat_marker@meta.data, Xexp < 10)

plot <- seurat_marker@meta.data %>% pivot_longer(c("Xexp", "CHR1exp", "CHR2exp", "UKexp"), 
                                         values_to = "percent_exp", names_to = "Chr_exp_type") %>% 
  as.data.frame %>% 
  ggplot(aes(x = percent_exp, colour = customclassif)) + geom_density() + 
  facet_grid(treatment ~ Chr_exp_type, scales ='free')

ggsave("plots/relative_exp_of_chrs_per_celltype.pdf", plot, width = 18)
#### DOING A PSEUDOBULK THING: probably not worth using or looking at lol ------


sce <- subset(seurat_marker, customclassif != 'Unknown') %>% 
  as.SingleCellExperiment(assay = "RNA")

agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("customclassif", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]

y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))
keep <- filterByExpr(y, group=agg_cell_types$treatment)
y <- y[keep,]
mds_data <-cpm(y, log = TRUE) %>% 
  plotMDS()


FindY1 <- y$counts %>% as.data.frame() %>% 
  mutate(gene = rownames(y$counts)) %>% 
  pivot_longer(colnames(y$counts), values_to = "reads", names_to = 'combsample') %>% 
  merge(y$samples[,c("customclassif", "sample")], by.x = 'combsample', by.y = 0) %>% 
  merge(coords_genes[,c("GENEID", "mid", "TXCHROM_name")], by.x = 'gene', by.y = 'GENEID')


FindY2 <- FindY1 %>% 
  group_by(gene, customclassif, mid, TXCHROM_name) %>% 
  summarise(mean_exp = mean(reads)) %>% 
  filter(TXCHROM_name != "unknown") %>% filter(customclassif == "Mature spermatids")

