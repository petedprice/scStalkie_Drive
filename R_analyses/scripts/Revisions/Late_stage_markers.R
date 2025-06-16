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

## MARKER GENES ##
markers <- read.xlsx("data/MANUSCRIPT/Drive Manuscript Supplementary Tables.xlsx", sheet = 2)[1:30,c(1:6)] 
row_order <- c("M", "CYSC", "EC", "LC", "C", "E", "G", "G, PS", "PS, SS", "PS, SS, ST", "PS, ST", "ST", "LST")
markers$Drosophila.cell.type <- factor(markers$Drosophila.cell.type, levels = row_order)
markers <- markers[order(markers$Drosophila.cell.type),]



new_markers <- c("Strip", "Dpy-30L2", "orb", "orb2", "Mst84D", "Pif", "Bug") %>% paste0(collapse = "|")
ls_markers <- ortholog_table %>% 
  filter(grepl(new_markers, consensus_gene)) %>% 
  mutate(REF_GENE_NAME = gsub("_", "-", REF_GENE_NAME))
ls_markers
DotPlot(seurat_final, features = c("g10249", unique(ls_markers$REF_GENE_NAME)))


seurat_final[,seurat_final$treatment == "SR"]  %>% 
  DotPlot(., features = c('gene-2553', "gene-4272", "gene-6897"))

novel_markers <- FindMarkers(seurat_final, 
            ident.1 = "Late spermatids", 
            ident.2 = "Early spermatids", 
            group.by = "celltype",
            logfc.threshold = 0,
            min.pct = 0,
            only.pos = T,
            assay = "RNA") %>% 
  #filter(p_val_adj < 1) %>% 
  mutate(gene = rownames(.))

View(novel_markers)
DotPlot(seurat_final, features = c("gene-4272", "gene-8980", "gene-6153", "STRG.7013"))

filter(ortholog_table, gsub("_", "-", REF_GENE_NAME) %in% novel_markers$gene) %>% 
  View()

DotPlot(seurat_final, novel_markers$gene, group.by = 'treatment')

