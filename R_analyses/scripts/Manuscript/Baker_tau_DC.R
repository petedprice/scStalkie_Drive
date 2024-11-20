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
tau_data <- openxlsx::read.xlsx('data/Baker2016_S2_Dataset_rev.xlsx')
Idents(seurat_final) <- "celltype"
levels(seurat_final) <- c("Muscle", "Early cyst", "Late cyst", 
                          "GSC/Spermatogonia", "Primary spermatocytes", 
                          "Secondary spermatocytes", "Spermatids")
seurat_final@meta.data <- seurat_final@meta.data %>%  
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids")))

DefaultAssay(seurat_final) <- "RNA"

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

tau <- ortholog_table %>% 
  dplyr::select(REF_GENE_NAME, consensus_gene) %>% 
  unique() %>% 
  merge(XA, by.x = 'REF_GENE_NAME', by.y = 'genes', all.y = T) %>% 
  merge(tau_data[,c('Symbol', 'Tau', 'Tissue.Specificity.I', 'Tissue.Specificity.II')], 
        by.x = 'consensus_gene', by.y = 'Symbol')


