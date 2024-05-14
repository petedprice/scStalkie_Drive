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


load("indata/RData/integrated_seurat.RData")
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
seurat_integrated$treatment = "ST" #which samples start with SR treatment is SR
seurat_integrated$treatment[grepl("sr", seurat_integrated$sample)] = "SR"


metadata <- seurat_integrated@meta.data
mutation_read <- function(file, metadata){
  samp <- str_split(file, '/', simplify = T) %>% c() %>% tail(.,1) %>% 
    str_split(., '_', simplify = T) %>% c() %>% head(.,1)
  mut_rate <- read.table(file, sep = " ", header = F)
  colnames(mut_rate) <- c("Barcode", "pos", "refn", "altn", "CHROM", "ref", "alt")
  mut_rate$depth <- mut_rate$refn + mut_rate$altn
  mut_rate$Barcode <- paste(samp, mut_rate$Barcode, sep = "_")
  #metadata_tmp <- filter(metadata, sample == samp & is.na(celltype) == F)
  metadata_tmp <- filter(metadata, sample == samp)
  
  mut_rate_sc <- mut_rate %>% 
    filter(depth > 0) %>% 
    merge(metadata[
      c('cells')],
      #is.na(metadata$sctype_labels) == F,c('cells', 'celltype')],
      by.x = 'Barcode', by.y = "cells")
  mut_rate_sc$genotype <- case_when(mut_rate_sc$refn == 0 & mut_rate_sc$altn > 1 ~"alt",
                                    mut_rate_sc$altn == 0 & mut_rate_sc$refn > 1~"ref",
                                    mut_rate_sc$refn > 1 & mut_rate_sc$altn > 1 ~ "het")
  mut_rate_sc$sample <- samp
  return(mut_rate_sc)
}
varfiles = list.files("indata/varcall/clus_filtered_var/", full.names = T, pattern = "snp_summarise.txt.gz")
mut_data <- lapply(varfiles, mutation_read, metadata = metadata) %>% 
  bind_rows()
metadata$cluster_no <- metadata$integrated_snn_res.0.1
hets <- mut_data  %>% 
  filter(is.na(genotype) == F) %>% 
  filter(CHROM != "Chr_X") %>% 
  merge(metadata, by.x = c('Barcode', 'sample'), by.y = c('cells', 'sample')) %>%
  group_by(cluster_no, sample, Barcode) %>% 
  filter(is.na(genotype) == F) %>%
  mutate(genotype = case_when(genotype %in% c("ref", "alt") ~ "hom",
                              genotype == "het" ~ "het")) %>% 
  summarise(heterozygosity = length(which(genotype == "het"))/length(genotype), 
            homozygosity = length(which(genotype == "hom"))/length(genotype), 
            n = n())


boxplots <- hets %>% 
  filter(n > 20) %>% 
  group_by(sample, cluster_no) %>% 
  summarise(heterozygosity = mean(heterozygosity), 
            n = n()) %>% 
  filter(n > 30) %>% 
  ggplot(aes(x = cluster_no, y = log(heterozygosity + 0.0001))) + geom_boxplot()

Idents(seurat_integrated) <- seurat_integrated$integrated_snn_res.0.1
DP <- DimPlot(seurat_integrated, label = T)

ggsave(boxplots, file = "plots/het_celltypes_boxplots.pdf", width = 15, height = 6)
ggsave(DP, file = "plots/het_celltypes_dimplot.pdf", width = 15, height = 8)





