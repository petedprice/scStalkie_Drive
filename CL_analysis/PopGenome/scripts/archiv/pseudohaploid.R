#!/usr/bin/env Rscript

install.packages("vcfR", lib = '.')
#install.packages("dplyr"), lib ='.')

library(vcfR, lib = '.')
library(memuse, lib = '.')
library(tidyverse)  #, lib = '.')

args = commandArgs(trailingOnly=TRUE)
vcf=args[1]
contig=args[2]
sample_info=read.csv(args[3])
out=args[4]

### GETTING SEX SAMPLES ###
vcf_whole <- read.vcfR(vcf)

Females <- intersect(sample_info$run_accession[sample_info$Sex == "F"], colnames(vcf_whole@gt))
Males <- c(intersect(sample_info$run_accession[sample_info$Sex == "M"], colnames(vcf_whole@gt)), "LTF")

### GETTING GT FUNCTION ###
get_gt <- function(samp, gt){
  gs <- str_split(gt[,samp], ":", simplify = T)[,1]
  df <- data.frame(genotype = gs)
  colnames(df) <- samp
  return(df)
}

### MALE GENOTYPES ###
just_gt_Males <- lapply(Males, get_gt, gt = vcf_whole@gt) %>% bind_cols() 

### FEMALE GENOTYPES ###
just_gt_Females <- lapply(Females, get_gt, gt = vcf_whole@gt) %>% bind_cols()

### FUNCTION TO CONVERT DIPLOID MALES AND AUTOSOMES TO PSEUDOHAPLOID ###
hap_the_dip <- function(samp, gt){
  gs <- str_split(gt[,samp], "/|\\|", simplify = T) %>% 
    as.data.frame()
  
  pseudo_dip <- lapply(gs, function(x)(return(paste0(x, "|", x)))) %>% 
    bind_cols()
  
  colnames(pseudo_dip) <- paste0(samp, "_hap", 1:2)
  return(pseudo_dip)
}

### GETTING THE PSEUDOHAPLOID MALE GT MATRIX ###
pseudo_haploid_males_gts <- lapply(Males, hap_the_dip, gt = just_gt_Males) %>% 
  bind_cols()

### GETTING THE PSEUDOHAPLOID FEMALE GT MATRIX ###
pseudo_haploid_female_gts <- lapply(Females, hap_the_dip, gt = just_gt_Females) %>%
      bind_cols()

### MERGING ###
pseudo_haploid_gts <- cbind(data.frame(FORMAT = rep("GT:AD:DP"))) %>% 
  cbind(pseudo_haploid_males_gts, pseudo_haploid_female_gts) %>% 
  as.matrix()
pseudo_haploid_gts[pseudo_haploid_gts == "0|0"] <- "0|0:10,0:10"
pseudo_haploid_gts[pseudo_haploid_gts == "1|1"] <- "1|1:0,10:10"
pseudo_haploid_gts[pseudo_haploid_gts == "2|2"] <- "2|2:0,10:10"
pseudo_haploid_gts[pseudo_haploid_gts == "3|3"] <- "3|3:0,10:10"
pseudo_haploid_gts[pseudo_haploid_gts == "4|4"] <- "4|4:0,10:10"

pseudo_haploid_gts[pseudo_haploid_gts == ".|."] <- ".|.:0,0:0"
pseudo_haploid_gts[pseudo_haploid_gts == "."] <- ".|.:0,0:0"
pseudo_haploid_gts[pseudo_haploid_gts == ".:0,0:0"] <- ".|.:0,0:0"



pseudohaploid_vcf <- vcf_whole
pseudohaploid_vcf@gt <- pseudo_haploid_gts

write.vcf(pseudohaploid_vcf, out)
