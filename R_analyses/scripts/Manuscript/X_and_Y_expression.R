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
Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 4)/
  colSums(seurat_final@assays$RNA$counts > 4)
nXgene <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 4)
nXprop_toX <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 4)/
  colSums(seurat_final@assays$RNA$counts[Xgenes,])

seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')
seurat_final <- AddMetaData(seurat_final, X_exp < 10, 'Y_cells')
seurat_final <- AddMetaData(seurat_final, nXgene, 'nXgene')
seurat_final <- AddMetaData(seurat_final, nXprop_toX, 'nXprop_toX')
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
FPnXprop_toX <- FeaturePlot(seurat_final, features = c('nXprop_toX'), pt.size = 0.1, 
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

BPFPnXprop_toX <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = nXprop_toX)) + geom_boxplot() +
  facet_wrap(~treatment) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

pdf("plots/X_exp_boxplots.pdf", height = 20, width = 20)
figure <- ggarrange(BPX_exp, BPY_cells, BPXprop, BPnXgene, BPnFeature, ncol = 1, nrow = 5, common.legend = T) 
annotate_figure(figure, bottom = "Cell Type") 
dev.off()

############## ############## ############## ############## ##############


pdf("plots/X_exp_densityplots.pdf", height = 20, width = 20)
seurat_final@meta.data %>% 
  ggplot(aes(x = nFeature_RNA, colour = treatment)) + geom_density() + 
  facet_wrap(~celltype) 
dev.off()

pdf("plots/Xprop_densityplots.pdf", height = 20, width = 20)
seurat_final@meta.data %>% 
  ggplot(aes(x = Xprop, colour = treatment)) + geom_density() + 
  facet_wrap(~celltype) 
dev.off()



################# MEETING WITH ALISON -------------

seurat_final@meta.data %>% 
  group_by(celltype, treatment) %>% 
  summarise(n = n(), 
            noX = sum(nXprop_toX < 0.2)) %>% 
  mutate(prop = noX/n)



#### X marker genes 
bunched  <- c("gene-9553")
CG8668 <- "gene-10423"
crc <- "PB.592"
mtgo <- "gene-10296"

X_markers <- c(bunched, CG8668, crc, mtgo)
DotPlot(seurat_final, features = X_markers, split.by = 'treatment')
FeaturePlot(seurat_final, features = X_markers, split.by = 'treatment')

CGexp <- seurat_final@assays$RNA$counts[CG8668,] > 0
crcexp <- seurat_final@assays$RNA$counts[crc,] > 0
mtgoexp <- seurat_final@assays$RNA$counts[mtgo,] > 0

seurat_final <- AddMetaData(seurat_final, CGexp, 'CG8668')
seurat_final <- AddMetaData(seurat_final, crcexp, 'crc')
seurat_final <- AddMetaData(seurat_final, mtgoexp, 'mtgo')


seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, fill = CG8668)) + geom_bar(position = 'dodge') + 
  facet_wrap(~treatment) 
seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, fill = crc)) + geom_bar(position = 'dodge') + 
  facet_wrap(~treatment) 
seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, fill = mtgo)) + 
  geom_bar(position = 'dodge') + 
  facet_wrap(~treatment)

FeaturePlot(seurat_final, features = c('CG8668', 'crc', 'mtgo'), 
            split.by = 'treatment')


###################################


seurat_final@meta.data %>% 
  #group_by(celltype, sample) %>% 
  ggplot(aes(x = celltype, y = Xprop, fill = treatment)) + geom_boxplot() +
  stat_compare_means( aes(label = ..p.signif..))
  #facet_wrap(~sample)
  



##########################
seurat_final@meta.data %>% 
  ggplot(aes(x = Xprop, colour = treatment)) + 
  geom_density() + 
  facet_wrap(~celltype)




