library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
rm(list = ls())

load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_no_cellcycle_markers.RData")
obj1 <- seurat_integrated

load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_newref.RData")
obj2 <- seurat_integrated

dim(obj1)
dim(obj2)

# cells <- colnames(obj1)[1:100]
# 
# obj1 <- subset(obj1, cells = cells)
# obj2 <- subset(obj2, cells = cells)

md1 <- obj1@meta.data[,c("cells", "mitoRatio", "nCount_RNA", "nFeature_RNA", "nFeature_SCT")]
md2 <- obj2@meta.data[,c("cells", "mitoRatio", "nCount_RNA", "nFeature_RNA", "nFeature_SCT")]

md <- merge(md1, md2, by = "cells")

(md$mitoRatio.x < md$mitoRatio.y) %>% sum()


missing_cells <- md2$cells[which(!md2$cells %in% md1$cells)]
length(missing_cells)

missing_md <- md2[md2$cells %in% missing_cells,]
View(missing_md)
View(md1)
View(md2)

