## based off http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
# load libraries and functions
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
library("biomaRt")
library(clusterProfiler)
library(ggrepel)
library(lme4)
library(GenomicFeatures)
library(ggrepel)
library(lme4)
library(AER)
library(nnet)
rm(list = ls())

load("data/RData/seurat_final.RData")
seurat_final$treatment <- "SR"
seurat_final$treatment[grep("st", seurat_final$sample)] <- "ST"
DefaultAssay(seurat_final) <- "integrated"
ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
celltype_table <- read.table("data/celltype_table.txt", header = T)


seurat_final <- subset(seurat_final, celltype != "Muscle")

germ <- c("GSC/Spermatogonia", "Primary spermatocytes", "Secondary spermatocytes", "Early spermatids", "Late spermatids")
cyst <- c("Early cyst", "Late cyst")
early_germ <- c("GSC/Spermatogonia", "Primary spermatocytes", "Secondary spermatocytes")
late_germ <- c("Late spermatids", "Early spermatids")

seurat_final@meta.data %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>%
  mutate(celltype2 = ifelse(celltype %in% germ, "Germ", "Cyst")) %>%
  mutate(binom = ifelse(celltype2 == "Germ", 0, 1)) %>%
  glmer(binom ~ treatment + (1|sample), family = "binomial", data = .) %>% 
  summary()

seurat_final@meta.data %>% 
  filter(celltype %in% germ) %>%
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>%
  mutate(celltype2 = ifelse(celltype %in% early_germ, "early_germ", "late_germ")) %>%
  mutate(binom = ifelse(celltype2 == "early_germ", 0, 1)) %>%
  glmer(binom ~ treatment + (1|sample), family = "binomial", data = .) %>% 
  summary()

seurat_final@meta.data %>% 
  filter(celltype %in% cyst) %>%
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>%
  mutate(binom = ifelse(celltype == "Early cyst", 0, 1)) %>%
  glmer(binom ~ treatment + (1|sample), family = "binomial", data = .) %>%
  summary()


