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

load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_newref_Seurat5.1.RData")
obj3 <- seurat_integrated


md1 <- obj1@meta.data[,c("cells", "mitoRatio", "nCount_RNA", "nFeature_RNA", "nFeature_SCT")]
md2 <- obj2@meta.data[,c("cells", "mitoRatio", "nCount_RNA", "nFeature_RNA", "nFeature_SCT")]
md3 <- obj3@meta.data[,c("cells", "mitoRatio", "nCount_RNA", "nFeature_RNA", "nFeature_SCT")]


md2_unique_cells <- md2$cells[which(!md2$cells %in% md3$cells)]
md3_unique_cells <- md3$cells[which(!md3$cells %in% md2$cells)]

md2_unique_cells[1]
md3_unique_cells[1]
View(md2)
View(md3)