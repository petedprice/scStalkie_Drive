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
library(lme4)
library(GenomicFeatures)
library(ggrepel)
library(lme4)


load("data/RData/seurat_final.RData")
DefaultAssay(seurat_final) <- "RNA"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  intersect(rownames(seurat_final))%>% unique()


##############################
X_exp <- PercentageFeatureSet(seurat_final, features = Xgenes)
A_exp <- PercentageFeatureSet(seurat_final, features = Agenes)
Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)/
  colSums(seurat_final@assays$RNA$counts > 0)
nXgene <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)


seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')
seurat_final <- AddMetaData(seurat_final, X_exp < 10, 'Y_cells')
seurat_final <- AddMetaData(seurat_final, nXgene, 'nXgene')

############## ############## Feature Plots ############## ############## 

UMAP <- DimPlot(seurat_final, label = T)
FPX_exp <- FeaturePlot(seurat_final, features = c('X_exp'), pt.size = 0.1, 
            split.by = 'treatment') + 
  theme(legend.position = c(0.1,0.2))
FPY_cells <- FeaturePlot(seurat_final, features = c('Y_cells'), pt.size = 0.1, 
            split.by = 'treatment') +
  theme(legend.position = c(0.1,0.2))
FPXprop <- FeaturePlot(seurat_final, features = c('Xprop'), pt.size = 0.1, 
            split.by = 'treatment') +
  theme(legend.position = c(0.1,0.2))
FPnXgene <- FeaturePlot(seurat_final, features = c('nXgene'), pt.size = 0.1, 
            split.by = 'treatment') +
  theme(legend.position = c(0.1,0.2))
FPnFeature <- FeaturePlot(seurat_final, features = c('nFeature_RNA'), pt.size = 0.1, 
            split.by = 'treatment') +
  theme(legend.position = c(0.1,0.2))

pdf("plots/X_exp_featureplots.pdf", height = 25, width = 20)
ggarrange(UMAP, FPX_exp, FPY_cells, FPXprop, FPnXgene, FPnFeature, ncol = 1, nrow = 6, 
          heights = c(2,1,1,1,1,1))
dev.off()

############## ############## ############## ############## ############## 




############## ############## BOX PLOTS ############## ############## 
BPX_exp <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = X_exp)) + geom_boxplot() +
  facet_wrap(~treatment) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

BPY_cells <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, fill = Y_cells)) + geom_bar(position = 'dodge') +
  facet_wrap(~treatment) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

BPXprop <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = Xprop)) + geom_boxplot() +
  facet_wrap(~treatment) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

BPnXgene <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = nXgene)) + geom_boxplot() +
  facet_wrap(~treatment) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

BPnFeature <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA)) + geom_boxplot() +
  facet_wrap(~treatment)

pdf("plots/X_exp_boxplots.pdf", height = 20, width = 20)
figure <- ggarrange(BPX_exp, BPY_cells, BPXprop, BPnXgene, BPnFeature, ncol = 1, nrow = 5, common.legend = T) 
annotate_figure(figure, bottom = "Cell Type") 
dev.off()

############## ############## ############## ############## ##############


pdf("plots/X_exp_densityplots.pdf", height = 20, width = 20)
seurat_final@meta.data %>% 
  ggplot(aes(x = X_exp, colour = treatment)) + geom_density() + 
  facet_wrap(~celltype) 
dev.off()

pdf("plots/Xprop_densityplots.pdf", height = 20, width = 20)
seurat_final@meta.data %>% 
  ggplot(aes(x = Xprop, colour = treatment)) + geom_density() + 
  facet_wrap(~celltype) 
dev.off()

