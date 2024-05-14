#!/usr/bin/env Rscript
##### LIBRARIES -------
library(Seurat)
library(tidyr)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(dplyr)
library(ggplot2)
#library(ggpubr)

###### SETTING UP INPUT COMMANDS ----
args = commandArgs(trailingOnly=TRUE)
outdatapath = paste(args[2], "/outdata/", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(args[2], "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)

##### FUNCTIONS ----
load(args[1])
samples=as.list(read.table(args[3]))
sample_info=read.table(args[4], header = F, sep = ",")
colnames(sample_info) <- c("sample", "species", 'ref', "sex", "stage")

mutation_read <- function(sample){
  samp = gsub("_", "", sample) %>% toupper()
  mut_rate <- read.table(paste(sample, "_mutation_data.txt.gz", sep = ""), sep = " ", header = F) 
  colnames(mut_rate) <- c("Barcode", "pos", "refn", "altn", "CHROM", "ref", "alt")
  mut_rate$Barcode <- paste(samp, mut_rate$Barcode, sep = "_")
  metadata_tmp <- filter(metadata, sample == samp & is.na(sctype_labels) == F)
  mut_rate_sc <- mut_rate %>% merge(metadata[
    is.na(metadata$sctype_labels) == F,c('cells', 'sctype_labels')], 
    by.x = 'Barcode', by.y = "cells")
  mut_rate_sc$sample <- sample
  return(mut_rate_sc)
}


metadata <- seurat_marker@meta.data


c=0
for (f in sample_info$sample){
  c=c+1
  print(f)
  if (c==1){
    merged_temp <- mutation_read(f)
  } else {
    temp <- mutation_read(f)
    merged_temp <- rbind(merged_temp, temp)
  }
}

print("done mutation_read")
merged_sc_mut <- merged_temp

somatic_snps <- filter(merged_sc_mut, sctype_labels != "germ") %>% 
  group_by(pos, CHROM, sample) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn)), 
          ref = ref[1], 
          alt = names(which.max(table(alt)))) %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref", 
                              refn > 0 & altn > 0 ~ "het")) 
print("done somatic snps")
#write.csv(somatic_snps,  file=gzfile("/fastdata/bop20pp/Avian_scRNAseq/R_analyses/mut_id/somatic_snps.csv.gz"))

somatic_snps_filt <- somatic_snps %>% 
  filter(depth > 4 & genotype != "het")
  
germ_snps <- filter(merged_sc_mut, sctype_labels == "germ") %>% 
  group_by(pos, CHROM, Barcode, sample, sctype_labels) %>% 
  reframe(refn = sum(refn), 
          altn = sum(altn), 
          depth = sum(c(altn, refn)), 
          ref = ref[1], 
          alt = names(which.max(table(alt)))) %>% 
  filter(pos %in% somatic_snps$pos) %>% 
  mutate(genotype = case_when(refn == 0 ~"alt", 
                              altn == 0 ~"ref",
                              (refn > 0 & altn > 0) ~"het"))
print("done germ snps")



#write.csv(germ_snps,  file=gzfile(
#  "/fastdata/bop20pp/Avian_scRNAseq/R_analyses/mut_id/germ_snps.csv.gz"))


germ_somatic <-  germ_snps %>% 
  merge(somatic_snps_filt, by = c('pos', 'sample', 'CHROM'), suffixes = c( '.other','.somatic')) %>% 
  mutate(mutated = case_when(genotype.other != genotype.somatic  & genotype.somatic != "het" ~ "mutation", 
                             genotype.other == genotype.somatic ~ "none"))

write.table(germ_somatic, paste(germ_somatic$sample[1], "_germ_somatic.txt", sep = ""))

print("doing gs file")


gs <- germ_somatic %>% 
  filter(genotype.other == 'het' | genotype.somatic != 'het') %>% 
  group_by(Barcode, sample, sctype_labels) %>% 
  reframe(nmuts = length(which(mutated == 'mutation')),
          nnone = length(which(mutated == "none")))
gs$mut_level <- gs$nmuts/(gs$nnone+gs$nmuts)


#write.table(gs, paste(gs$sample[1], "_gs.txt", sep = ""))

#seurat_mutant <- seurat_marker
mutant_metadata <- merge(metadata, gs[,c(1,4,5,6)], by.x = "cells", by.y = "Barcode", all.x = F)

#save(seurat_mutant, file = paste(outdatapath, "/seurat_mutant.RData", sep = ""))
write.csv(mutant_metadata, paste(mutant_metadata$sample[1], "mutant_metadata.csv", sep = "_"))
