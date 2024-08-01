library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(gridExtra)
library(data.table)
rm(list = ls())

load("data/RData/seurat_final.RData")
load("data/RData/tau_data.RData")
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

samples <- c("sr1", "sr2", "sr3", "sr5", "st1","st2","st3","st5")

run_checks = "yes"
if (run_checks == "yes"){
  
  metadata <- seurat_final@meta.data

  
  #thresholds <- data.frame(het = c(0,1,2,3), hom = c(0,2,4,6), nsnps = rep(c(0,5,10), each = 4))
  thresholds <- data.frame(het = c(0), hom = c(1), nsnps = rep(c(5), each = 1))
  metadatas <- list()
  het_data_sum_save <- list()
  het_data_save <- list()
  c=0
  ot <- ortholog_table[,c("chr", "REF_GENE_NAME", "start", "end")] %>% 
    rename(Chr = chr) %>% 
    #filter(Chr != "Chr_X") %>% 
    filter(REF_GENE_NAME %in% tau_genes$genes) %>% 
    merge(tau_genes[,c("genes", "tau_type")], by.x = "REF_GENE_NAME", by.y = "genes") %>% 
    as.data.table()
  samp = samples[1]
  for (samp in samples){  
    print(samp)
    snps <- read.table(paste0("data/scAlleleCount/stringent/", samp, "_snps.txt")) %>% 
      #filter(V1 != "Chr_X") %>% 
      mutate(Pos2 = V2) %>% as.data.table()
      
    snps$snp_indx <- 1:nrow(snps)
    colnames(snps) <- c("Chr", "Pos", "Ref", "Alt", "Pos2", "snp_indx")
    setkey(snps, Chr, Pos, Pos2)
    setkey(ot, Chr, start, end)
    #ot[, c("start", "end") := .(start, end)]
    snps_tau <- foverlaps(snps, ot, by.x = c("Chr", "Pos", "Pos2"), 
                        by.y = c("Chr", "start", "end"), type = "within", nomatch = 0)
    snps_tau <- snps_tau[,c("Chr", "Pos", "Ref", "Alt", "snp_indx", "REF_GENE_NAME", "tau_type")]
    
    
    barcodes <- read.table(paste0("data/scAlleleCount/stringent//data/", samp, "_barcodes.tsv"))[,1]
    cov <- readMM(file = paste0("data/scAlleleCount/stringent/data/", samp, "covmat.mtx"))
    alt <- readMM(file = paste0("data/scAlleleCount/stringent/data/", samp, "altmat.mtx"))
    ref <- readMM(file = paste0("data/scAlleleCount/stringent/data/", samp, "refmat.mtx"))
    
    #subsetting matrices for cells we have present in our seurat object
    cells <- match(paste0(samp, "_", barcodes), rownames(seurat_final@meta.data))
    cell_names <- rownames(seurat_final@meta.data)[cells]
    cell_names <- cell_names[!is.na(cell_names)]
    cn_df <- data.frame(cell_barcode = cell_names, cell = 1:length(cell_names))
    cell_indx <- which(!is.na(cells))
    covsub <- cov[cell_indx,]
    altsub <- alt[cell_indx,]
    refsub <- ref[cell_indx,]
    
    #summarising the matrices 
    scov <- summary(covsub) %>% 
      filter(j %in% snps_tau$snp_indx)
    salt <- summary(altsub) %>% 
      filter(j %in% snps_tau$snp_indx)
    sref <- summary(refsub) %>% 
      filter(j %in% snps_tau$snp_indx)
    
    #merging with i bring the cell and j being the position and x and y being the alt number and cov being total number (coverage)
    het_data <- merge(salt, scov, by=c("i", "j"), all = T) %>% 
      # change NA values to 0
      replace_na(list(x.x = 0, x.y = 0)) %>% 
      dplyr::rename(cell = i, snp_indx = j, altn = x.x, covn = x.y) %>%
      merge(snps_tau, by = 'snp_indx') %>% 
      mutate(homozygosity = altn/covn) %>% 
      merge(cn_df, by = 'cell') %>% #add column which is cov -alt 
      mutate(refn = covn - altn)
    het_data$sample <- samp
    het_data_save[[samp]] <- het_data
    for (tt in 1:nrow(thresholds)){
      print(tt)
      c=c+1
      het_data_sum <- het_data %>% 
        group_by(cell_barcode, Chr, tau_type) %>% 
        summarise( #summarise so hom is no of pos where homozygosity is 1 and het is no of pos where homozygosity is not 1
          homn = length(which(homozygosity %in% c(0,1) & covn > thresholds[tt,2])), 
          homd= sum(covn[homozygosity %in% c(0,1) & covn > thresholds[tt,2]]),
          hetn = length(which(homozygosity > 0 & homozygosity < 1 & altn > thresholds[tt,1] & refn > thresholds[tt,1] & covn > thresholds[tt,2])),
          hetd= sum(covn[!homozygosity %in% c(0,1) & covn > thresholds[tt,2]]),
          total_hethoms = (hetn + homn), 
          hethomd=homd+hetd,
          readsn = sum(covn),
          snpsn = length(unique(snp_indx)),
          snpcov = sum(covn)/snpsn) %>% 
        mutate(heterozygosity = homn/total_hethoms)
      het_data_sum$coverage_check <- "OK"
      het_data_sum$coverage_check[het_data_sum$total_hethoms < thresholds[tt,3] | is.na(het_data_sum$snpsn)] <- "low_coverage"
      het_data_sum$sample <- samp
      het_data_sum$depth_threshold <- thresholds[tt,1]
      het_data_sum$snp_threshold <- thresholds[tt,3]
      
      
      het_data_sum_save[[c]] <- het_data_sum
    }
  }
  
  
  
  new_meta_data <- bind_rows(het_data_sum_save) %>% 
    merge(metadata, by.x = c('cell_barcode', 'sample'), by.y = c('cells', 'sample'), all.y = F)
  new_meta_data$ploidy_class <- ifelse(new_meta_data$heterozygosity < 0.80, "diploid", "haploid")
  new_meta_data$ploidy_class[is.na(new_meta_data$ploidy_class)] <- "Missing"
  
  summary_table <- new_meta_data %>% 
    group_by(depth_threshold, snp_threshold, celltype, ploidy_class, tau_type, Chr) %>%
    summarise(n =n())
  grid.table(summary_table)
  
  #save(summary_table, het_data_save, new_meta_data, file = "data/scAlleleCount//ploidy_params_check.RData")
} else if (run_checks == "no"){
  load(file = "data/scAlleleCount//ploidy_params_check.RData")
}



bar_graph <- summary_table %>% 
  ggplot(aes(x = celltype, fill = ploidy_class, y = n)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~tau_type + Chr, nrow = 4) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
wbar_graph


ggsave("plots/ploidy_check.png", bar_graph, width = 40, height = 20, units = "in")

write.table(summary_table, 'data/scAlleleCount/ploidy_check_summary.tsv', sep = '\t', quote = F, row.names = F)
####


