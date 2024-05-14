## based off http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
# load libraries and functions
library(dplyr)
library(Seurat)
library(openxlsx)
library(ggpubr)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(lme4)
library(vcfR)


load("indata/RData/integrated_seurat.RData")
vcfs <- list.files("indata/varcall/filtered_vcfs_clusters/", full.names = T, pattern = "vcf")

vcf_read <- function(x){
  vcf <- read.vcfR(x, verbose = FALSE)
  if (nrow(vcf@gt) > 0) {
    out <- cbind(vcf@fix, vcf@gt) %>% as.data.frame()
    out$sample <- basename(x) %>% 
      str_split("_", simplify = T) %>% .[,1]
    out$barcode <- basename(x) %>% 
      str_split("_", simplify = T) %>% .[,2] %>% 
      gsub("clus", "", .)
    out$genotype <- str_split(out[,10], ":", simplify = T) %>% .[,1]
    colnames(out)[10] <- "OUT"
    dpi <- lapply(out[,9], function(x)(which(str_split(x, ":", simplify = T) == "DP"))) %>% unlist()
    adi <- lapply(out[,9], function(x)(which(str_split(x, ":", simplify = T) == "AD"))) %>% unlist()
    
    out$depth <- lapply(1:nrow(out), function(x)(str_split(out[x,10], ":", simplify = T)[dpi[x]])) %>% unlist() %>% as.numeric()
    out$allele_depth <- lapply(1:nrow(out), function(x)(str_split(out[x,10], ":", simplify = T)[adi[x]])) %>% unlist() 
    
    min_ad <- function(x){
      dps <- as.numeric(str_split(x, ",", simplify = T))
      return(min(dps[dps > 0]))
    }
    out$min_ad <- lapply(out$allele_depth, min_ad) %>% unlist()
    out <- out %>% filter(genotype != "./.") %>% 
      filter(min_ad > 1)

    if (nrow(out) > 0){
      return(out)
    }
  }
}


vcf <- sapply(vcfs, vcf_read) %>% bind_rows()


vcf_summary <- vcf %>% 
  group_by(barcode, CHROM, sample, genotype) %>% #get counts for each type of genotype 
  filter(genotype != "./.") %>%
  summarise(count = n())

vcf_summary <- vcf %>% 
  group_by(barcode, CHROM, sample, genotype) %>% #get counts for each type of genotype 
  #filter(CHROM != "Chr_X") %>% 
  group_by(barcode, sample) %>% 
  summarise(heterozygosity = sum(genotype %in% c("0/1", "1/2"))/sum(genotype != "./."),
            homozygosity = sum(genotype == "1/1")/sum(genotype != "./."), 
            n = n())


vcf_summary_md <- vcf_summary %>% 
  filter(n > 2) %>% 
  mutate(samp_barcode = paste0(sample, "_", barcode)) %>%
  merge(seurat_integrated@meta.data[,c("cells", "integrated_snn_res.0.1")], by.x= "samp_barcode", by.y = "cells")


boxp <- vcf_summary_md %>% 
  ggplot(aes(x = integrated_snn_res.0.1, y = heterozygosity, fill = integrated_snn_res.0.1)) + 
  geom_boxplot()
dimp <- DimPlot(seurat_integrated, group.by = "integrated_snn_res.0.1", label = T)

ggarrange(boxp, dimp, ncol = 2)
table(vcf_summary_md$integrated_snn_res.0.1)
