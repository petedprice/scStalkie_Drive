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
ortholog_table$length <- ortholog_table$end - ortholog_table$start
ortholog_table$mid <- ortholog_table$start + (ortholog_table$length/2)


pg_markers <- readxl::read_excel("data/PopGroup_markers.xlsx") 
colnames(pg_markers)[1] <- "seurat_markers"

for (i in 1:nrow(pg_markers)){
  seurat_marker$celltype[seurat_marker$seurat_clusters %in% pg_markers$seurat_markers[i]] <- pg_markers$New_consensus[i]
}
seurat_marker <- subset(seurat_marker, celltype != 'Unknown' & celltype != "Accessory Gland?" & celltype != "Muscle")
seurat_marker$celltype[seurat_marker$`integrated_snn_res.1.5` == 32] <- "Spermatids"


############### VARIATION DATA -------------
metadata <- seurat_integrated@meta.data
file = "indata/varcall/filtered_var300124/filtered_var/st1_Chr_X_snp_summarise.txt.gz"
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
  mut_rate_sc$genotype <- case_when(mut_rate_sc$refn == 0 ~"alt",
                                    mut_rate_sc$altn == 0 ~"ref",
                                    mut_rate_sc$refn > 0 & mut_rate_sc$altn > 0 ~ "het")
  mut_rate_sc$sample <- samp
  return(mut_rate_sc)
}

varfiles <- list.files("indata/varcall/filtered_var190124/", full.names = T)
varfiles <- varfiles[grep("Chr_X", varfiles)]
st1_X_mutdata <- mutation_read(
  "indata/varcall/filtered_var140224//st1_snp_summarise.txt.gz", 
  metadata) 
write.table(st1_X_mutdata, "indata/varcall/st1_compiled.csv", sep = ",", quote = F, row.names = F)
mut_data <- lapply(varfiles, mutation_read, metadata = metadata) %>% 
  bind_rows()
hd <- mut_rate_sc[which(mut_rate_sc$depth == max(mut_rate_sc$depth)),]$pos
hd_gene <- ortholog_table[which(ortholog_table$start < hd & ortholog_table$end > hd),]
seurat_marker@assays$RNA$counts["PB.7527",'sr1_GAAGCGAGTTGGAGAC-1']
summaryX <- mut_data %>% filter(CHROM == "Chr_X") %>% 
  filter(celltype %in% c("Epithelial", "Cyst", "Muscle cells")) %>% 
  group_by(pos, sample) %>% 
  summarise(refn = sum(refn), altn = sum(altn)) %>%
  mutate(genotype = case_when(refn == 0 ~"alt",
                              altn == 0 ~"ref",
                              refn > 5 & altn > 5 ~ "het"))

summaryX %>% 
  filter(is.na(genotype) ==F) %>% 
  ggplot(aes(x = pos, colour = genotype)) + geom_density()

summaryX %>% 
  filter(is.na(genotype) ==F) %>% 
  ggplot(aes(x = pos, y = genotype)) + geom_point(alpha = 0.01)

mut_data$celltype2 <- NA
mut_data$celltype2[mut_data$celltype %in% c("Cyst", "Muscle cells")] <- "Somatic"
mut_data$celltype2[mut_data$celltype %in% c("Spermatids")] <- "Germ"

metadata_clus <- metadata %>% 
  dplyr::select(`integrated_snn_res.0.1`, cells) %>% 
  rename(cluster_no =`integrated_snn_res.0.1`)

hets <- st1_X_mutdata %>% #filter(celltype2 %in% c("Somatic", "Germ")) %>% 
  filter(is.na(genotype) == F) %>% 
  merge(metadata_clus, by.x = 'Barcode', by.y = 'cells') %>%
  group_by(cluster_no, sample, CHROM) %>% 
  filter(is.na(genotype) == F) %>%
  mutate(genotype = case_when(genotype %in% c("ref", "alt") ~ "hom",
                              genotype == "het" ~ "het")) %>% 
  summarise(heterozygosity = length(which(genotype == "het"))/length(genotype), 
            homozygosity = length(which(genotype == "hom"))/length(genotype), 
            n = n())


clus_order <- c(12,5,30,0,6,1,10,14,24,29,23,26,27,
                19,2,3,18,8,17,16,20,7,11,13,28,21,25,4,9,15,22)
hets$cluster_no <- factor(hets$cluster_no, levels = clus_order)
boxplots <- hets %>% 
  ggplot(aes(x = cluster_no, y = heterozygosity, fill = factor(cluster_no))) + geom_boxplot() +
  facet_wrap(~CHROM, ncol = 2) + theme_bw()

boxplots <- hets %>% 
  ggplot(aes(x = cluster_no, y = heterozygosity, fill = factor(CHROM))) + geom_boxplot() +
   theme_bw()

DP <- DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.1', label = T)

ggsave(boxplots, file = "plots/het_celltypes_boxplots.pdf", width = 15, height = 6)
ggsave(DP, file = "plots/het_celltypes_dimplot.pdf", width = 15, height = 8)

hets %>% 
  ggplot(aes(x = CHROM, y = heterozygosity, fill = celltype2)) + 
  geom_bar(stat = "identity", position = "dodge")

metadatahets <- merge(metadata, hets, by.x = c('cells', "sample", "celltype"), 
                      by.y = c('Barcode', "sample", "celltype"), all.x = T)
metadatahets$heterozygosity[is.na(metadatahets$heterozygosity)] <- 0.00001
seurat_marker <- AddMetaData(seurat_marker, metadatahets)
seurat_marker@meta.data <- metadatahets
FeaturePlot(seurat_marker, features = c('heterozygosity'), split.by = 'sample')


metadatahets %>% 
  ggplot(aes(x = celltype, y = heterozygosity * 100, fill = treatment)) + geom_boxplot()


metadatahets %>% 
  ggplot(aes(x = celltype, y = log10GenesPerUMI, fill = treatment)) + geom_boxplot()





################## vcfs include 
sp = 'st1'
mut_data
keep_snps <- function(sample){
  snps <- filter(sample_summaries, sample == sp) %>% 
    dplyr::select(POS) %>% c() %>% unlist()
  md_ss <- filter(mut_data, sample == sp & pos %in% snps)
  return(md_ss)
}

filt2 <- lapply(unique(mut_data$sample), keep_snps) %>% bind_rows()


filt2 %>% 
  group_by(sample, CHROM, celltype) %>%
  summarise(heterozygosity = length(which(genotype == "het"))/length(genotype), 
            homozygosity = length(which(genotype %in% c("ref", "alt")))/length(genotype), 
            n = n()) %>%
  ggplot(aes(x = CHROM, y = heterozygosity, fill = celltype, colour = sample)) +
  geom_bar(stat = "identity", position = "dodge")
  




