################# load libraries and functions and data prep ###################
library(dplyr)
library(Seurat)
library(openxlsx)
library(ggpubr)
library(scuttle)
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(ggrepel)
library(GenomicFeatures)
library(ggrepel)
library(cowplot)
rm(list = ls())

#### LOAD DATA AND PREP DATA ---
load("data/RData/seurat_final.RData")
load("data/RData/DEG_DC.RData")
load("data/trajectory/sce_GAMed_8k.RData")
load("data/RData//traj_analyis.RData")

seurat_final <- JoinLayers(seurat_final)
Idents(seurat_final) <- "celltype"
levels(seurat_final) <- c("Muscle", "Early cyst", "Late cyst", 
                          "GSC/Spermatogonia", "Primary spermatocytes", 
                          "Secondary spermatocytes", "Early spermatids", "Late spermatids")
seurat_final@meta.data <- seurat_final@meta.data %>%  
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids")))

DefaultAssay(seurat_final) <- "RNA"

ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")


ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]




Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()

nFeatures_RNA_autosomes <- colSums(seurat_final@assays$RNA$counts[Agenes,] > 0)

seurat_final <- AddMetaData(seurat_final, nFeatures_RNA_autosomes, 'nFeature_RNA_autosomes')

X_exp <- PercentageFeatureSet(seurat_final, features = Xgenes)
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')

Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)/
  colSums(seurat_final@assays$RNA$counts > 1)
Xgene_prop2 <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)/
  nrow(ortholog_table[ortholog_table$chr == "Chr_X",])



seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')
seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop2, col.name = 'Xprop2')


seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = Xprop, fill = celltype)) + geom_boxplot()

