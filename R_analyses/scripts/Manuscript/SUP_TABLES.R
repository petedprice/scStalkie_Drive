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



load("data/RData/seurat_final.RData")
DefaultAssay(seurat_final) <- "RNA"
Idents(seurat_final) <- "celltype"
seurat_final <- JoinLayers(seurat_final)
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

### MARKERS using FindAllMarkers ###
markers <- seurat_final %>% 
  subset(., treatment == "ST")  %>%
  FindAllMarkers(., only.pos = T, min.pct = 0.25, min.logfc = 2)


markers %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.01 & pct.2 < 0.1) %>% 
  mutate(cluster = factor(cluster, levels = c("Muscle", "Early cyst",
                                              "Late cyst", "GSC/Spermatogonia", 
                                              "Primary spermatocytes", "Secondary spermatocytes", 
                                              "Early spermatids", "Late spermatids"))) %>%
  group_by(cluster) %>% 
  top_n(., 20, -p_val_adj) %>% 
  merge(unique(ortholog_table[,c("REF_GENE_NAME", "consensus_gene")]), by.x = 'gene', by.y = 'REF_GENE_NAME') %>% 
  mutate(consensus_gene = ifelse((consensus_gene == gene), NA, consensus_gene)) %>% 
  arrange(cluster, p_val_adj) %>% 
  mutate(pct.1 = 100 * pct.1, 
         pct.2 = 100 * pct.2) %>%
  rename(
    `Reference Gene` = gene, 
    log2FC = avg_log2FC,
    `% of cells expressing in ref cluster` = pct.1,
    `% of cells expressing in other clusters` = pct.2,
    Cluster = cluster,
    `Drosophila ortholog` = consensus_gene
  ) %>% 
  write.csv("data/MANUSCRIPT/FindAllMarkers_ST.csv", row.names = F)


seurat_final@meta.data %>% 
  group_by(sample, treatment) %>% 
  rename(Sample = sample, 
         type = treatment, 
         `Cell type` = celltype) %>%
  summarise(`Number of cells` = n(), 
            `Median no of features` = mean(nFeature_RNA), 
            `Median no of UMIs` = (median(nCount_RNA))) %>% 
  write.csv("data/MANUSCRIPT/Cell_numbers_samptreat.csv", row.names = F)

seurat_final@meta.data %>% 
  group_by(sample, treatment, celltype) %>% 
  rename(Sample = sample, 
         type = treatment, 
         `Cell type` = celltype) %>%
    summarise(`Number of cells` = n(), 
            `Median no of features` = median(nFeature_RNA), 
            `Median no of UMIs` = median(nCount_RNA)) %>% 
  write.csv("data/MANUSCRIPT/Cell_numbers_samptreatcell.csv", row.names = F)
 


