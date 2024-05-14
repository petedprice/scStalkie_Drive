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
seurat_marker$celltype = seurat_marker$sctype_labels
seurat_marker$treatment = "ST" #which samples start with SR treatment is SR
seurat_marker$treatment[grepl("sr", seurat_marker$sample)] = "SR"


############### VARIATION DATA -------------
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
varfiles = list.files("indata/varcall/filtered_var/", full.names = T, pattern = "snp_summarise.txt.gz")
mut_data <- lapply(varfiles, mutation_read, metadata = metadata) %>% 
  bind_rows()
metadata$cluster_no <- metadata$integrated_snn_res.0.2
hets <- mut_data %>% #filter(celltype2 %in% c("Somatic", "Germ")) %>% 
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







max_het_bc <- hets$Barcode[which.max(hets$heterozygosity)]
filter(mut_data, Barcode == max_het_bc) %>% 
  filter(!is.na(genotype) & CHROM != "Chr_X") %>% 
  group_by(Barcode) %>% 
  filter(is.na(genotype) == F) %>%
  mutate(genotype = case_when(genotype %in% c("ref", "alt") ~ "hom",
                              genotype == "het" ~ "het")) %>% 
  summarise(heterozygosity = length(which(genotype == "het"))/length(genotype), 
            homozygosity = length(which(genotype == "hom"))/length(genotype), 
            n = n())
  
filter(hets, Barcode == max_het_bc)




boxplots <- hets %>% 
  filter(n > 20) %>% 
  ggplot(aes(x = cluster_no, y = log(heterozygosity + 0.0001))) + geom_boxplot()

Idents(seurat_integrated) <- seurat_integrated$integrated_snn_res.0.2
DP <- DimPlot(seurat_integrated, label = T)

ggsave(boxplots, file = "plots/het_celltypes_boxplots.pdf", width = 15, height = 6)
ggsave(DP, file = "plots/het_celltypes_dimplot.pdf", width = 15, height = 8)

hets %>% 
  ggplot(aes(x = CHROM, y = heterozygosity, fill = celltype2)) + 
  geom_bar(stat = "identity", position = "dodge")

metadatahets <- merge(metadata, hets, by.x = c('cells', "sample"), 
                      by.y = c('Barcode', "sample"), all.x = T)
metadatahets$heterozygosity[is.na(metadatahets$heterozygosity)] <- 0.00001
seurat_integrated <- AddMetaData(seurat_integrated, metadatahets)
seurat_integrated@meta.data <- metadatahets
FeaturePlot(seurat_integrated, features = c('heterozygosity'), split.by = 'sample')


metadatahets %>% 
  filter(!is.na(homozygosity)) %>% 
  filter(n > 5) %>% 
  ggplot(aes(x = cluster_no, y = homozygosity, fill = factor(cluster_no))) + geom_boxplot()
  
  
hets %>% 
  filter(!is.na(homozygosity)) %>% 
  filter(n > 50) %>% 
  group_by(cluster_no) %>%
  summarise(mean_het = mean(heterozygosity) * 100,
            var = var(heterozygosity),
            n = n())




spermatid_barcodes <- filter(metadata, 
                             integrated_snn_res.0.2 == 0 & 
                               sample == 'st1') %>% 
  select(cells) %>% 
  rownames() %>% 
  gsub("st1_", "", .)
write.table(spermatid_barcodes, file = "data/st1_spermatid_barcodes.txt", quote = F, row.names = F, col.names = F)

muscle_barcodes <- filter(metadata, 
                             integrated_snn_res.0.2 == 15 & 
                               sample == 'st1') %>% 
  select(cells) %>% 
  rownames() %>% 
  gsub("st1_", "", .)
write.table(muscle_barcodes, file = "data/st1_muscle_barcodes.txt", quote = F, row.names = F, col.names = F)


