## based off http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(SCINA)
library(ggpubr)
#library(scater)
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

### FUNCTIONS ---
mutation_read <- function(file, metadata){
  samp <- str_split(file, '/', simplify = T) %>% c() %>% tail(.,1) %>% 
    str_split(., '_', simplify = T) %>% c() %>% head(.,1)
  mut_rate <- read.table(file, sep = " ", header = F)
  colnames(mut_rate) <- c("Barcode", "pos", "refn", "altn", "CHROM", "ref", "alt")
  mut_rate$depth <- mut_rate$refn + mut_rate$altn
  mut_rate$Barcode <- paste(samp, mut_rate$Barcode, sep = "_")
  metadata_tmp <- filter(metadata, sample == samp & is.na(celltype) == F)
  mut_rate_sc <- mut_rate %>% 
    filter(depth > 0) %>% 
    merge(metadata[
      is.na(metadata$sctype_labels) == F,c('cells', 'celltype')],
      by.x = 'Barcode', by.y = "cells")
  mut_rate_sc$genotype <- case_when(mut_rate_sc$refn == 0 ~"alt",
                                    mut_rate_sc$altn == 0 ~"ref",
                                    mut_rate_sc$refn > 0 & mut_rate_sc$altn > 0 ~ "het")
  mut_rate_sc$sample <- samp
  return(mut_rate_sc)
}


add_gene_func <- function(pos, ot = ortholog_table){
  gene <- ot[ot[,6] <= pos & ot[,7] >= pos,1]
  return(gene[1])
}

########## WORK 

########

#########################


######### DATA ------
load("data/RData/marker_seurat.RData")
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table$REF_GENE_NAME <- gsub("gene_", "gene-", ortholog_table$REF_GENE_NAME)
seurat_marker$celltype = seurat_marker$sctype_labels
seurat_marker$treatment = "ST" #which samples start with SR treatment is SR
seurat_marker$treatment[grepl("sr", seurat_marker$sample)] = "SR"
ortholog_table$length <- ortholog_table$end - ortholog_table$start
ortholog_table$mid <- ortholog_table$start + (ortholog_table$length/2)


pg_markers <- readxl::read_excel("data/PopGroup_markers.xlsx") 
colnames(pg_markers)[1] <- "seurat_markers"

for (i in 1:nrow(pg_markers)){
  seurat_marker$celltype[seurat_marker$seurat_clusters %in% pg_markers$seurat_markers[i]] <- pg_markers$New_consensus[i]
}
metadata <- seurat_marker@meta.data
varfiles <- list.files("indata/varcall/filtered_var190124/filtered_var", full.names = T)
varfiles <- varfiles[grep("Chr_X", varfiles)]

sr1_X_old <- read.table(
  "indata/varcall/filtered_var190124/filtered_var/sr1_Chr_X_snp_summarise.txt.gz") 
colnames(sr1_X_old) <- c("Barcode", "pos", "refn", "altn", "CHROM", "ref", "alt")

sr1_X_new <- read.table("indata/varcall/filtered_var300124/sr1_snp_summarise.txt.gz")
colnames(sr1_X_new) <- c("Barcode", "pos", "refn", "altn", "CHROM", "ref", "alt")

### ORGANISING DATA ----
sr1o <- sr1_X_old %>% 
  filter(CHROM == "Chr_X") %>% 
  group_by(pos, CHROM) %>% 
  summarise(depth = sum(refn + altn)) %>% 
  apply(., 1, add_gene_func)
geneo <- lapply(sr1o$pos, add_gene_func)
sr1o$gene <- unlist(geneo)
sr1of <- sr1o %>% group_by(gene) %>% 
  filter(is.na(gene) == F) %>% 
  summarise(md = mean(depth)) %>% as.data.frame()


sr1n <- sr1_X_new %>% 
  filter(CHROM == "Chr_X") %>% 
  group_by(pos, CHROM) %>% 
  summarise(depth = sum(refn + altn))
genen <- lapply(sr1n$pos, add_gene_func)
sr1n$gene <- unlist(genen)
sr1nf <- sr1n %>% group_by(gene) %>% 
  filter(is.na(gene) == F) %>% 
  summarise(md = mean(depth)) %>% as.data.frame()

sr1 <- subset(seurat_marker, sample == 'sr1')
mean_exp <- sr1@assays$RNA@counts %>% rowMeans() %>% 
  as.data.frame()
mean_exp$gene <- rownames(mean_exp)
colnames(mean_exp) <- c("mean_exp", "gene")
rownames(mean_exp) <- NULL
mean_exp <- as.data.frame(mean_exp)
####################

### DATA MERGING ----

merged_data <- sr1nf %>% 
  left_join(sr1of, by = "gene", suffix = c("_new", "_old"))
merged_data <- merged_data %>% 
  left_join(mean_exp, by = "gene")

plots <- list()
plots[[1]] <- merged_data %>% 
  ggplot(aes(x = log(md_old), y = log(mean_exp))) + geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm") + stat_cor()
plots[[2]] <- merged_data %>%
  ggplot(aes(x = log(md_new), y = log(mean_exp))) + geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm") + stat_cor()

plots[[3]] <- merged_data %>% 
  ggplot(aes(x = log(md_new), y = log(md_old))) + geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm") + stat_cor()

plots
### TO DO 
#READ IN OLD VARCALL DATA FOR SAMPLE SR1 and sumarise by position/chr 
#READ IN NEW VARCALL DATA FOR SAMPLE SR1 and sumarise by position/chr
#CALCULATE MEAN EXPRESSION FOR EACH GENE FOR SR1 for seurat_marker





