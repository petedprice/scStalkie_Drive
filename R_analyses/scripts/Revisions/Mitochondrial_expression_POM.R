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
load("data/RData/DEG_DC.RData")

seurat_final <- JoinLayers(seurat_final)

seurat_final$celltype[seurat_final$celltype %in% c("Early spermatids", "Late spermatids")] <- "Spermatids"
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
seurat_final$treatment <- "SR"
seurat_final$treatment[grep("st", seurat_final$sample)] <- "ST"
ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene, FBconcensus) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)


POM_MITO <- seurat_final@meta.data %>% 
  mutate(celltype = factor(celltype, levels <- c("Muscle", "Early cyst", "Late cyst", 
                            "GSC/Spermatogonia", "Primary spermatocytes", 
                            "Secondary spermatocytes", "Early spermatids", "Late spermatids"))) %>% 
  ggplot(aes(x = celltype, y = 100*mitoRatio, fill = treatment)) + 
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + 
  labs(fill = "", x = "", y = "% of transcripts mapping from mitochondrial genes") +
  scale_fill_grey(start = 0.35, end = 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black"))
ggsave("plots/POM_MITO.pdf", height = 6, width = 6)
system("open plots/POM_MITO.pdf")

