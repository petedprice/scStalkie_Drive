library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(gridExtra)
library(cowplot)
rm(list = ls())
load("data/RData/seurat_final.RData")
samples <- c("sr1", "sr2", "sr3", "sr5", "st1","st2","st3","st5")



  
metadata <- seurat_final@meta.data
mut_data <- list()
for (samp in unique(metadata$sample)){
  print(samp)

  
  snps <- read.table(paste0("data/scAlleleCount/relaxed/send/", samp, "_filt.recode.bed"))
  snps$snp_indx <- 1:nrow(snps)
  colnames(snps) <- c("Chr", "Pos", "Ref", "Alt", "snp_indx")
  barcodes <- read.table(paste0("data/scAlleleCount/relaxed//send/", samp, "_barcodes.tsv"))[, 1]
  cov <- readMM(file = paste0("data/scAlleleCount/relaxed//send/", samp, "covmat.mtx"))
  alt <- readMM(file = paste0("data/scAlleleCount/relaxed/send/", samp, "altmat.mtx"))
  ref <- readMM(file = paste0("data/scAlleleCount/relaxed/send/", samp, "refmat.mtx"))
  
  #subsetting matrices for cells we have present in our seurat object
  cells <- match(paste0(samp, "_", barcodes), rownames(seurat_final@meta.data))
  cell_names <- rownames(seurat_final@meta.data)[cells]
  cell_names <- cell_names[!is.na(cell_names)]
  cn_df <- data.frame(cell_barcode = cell_names, cell = 1:length(cell_names))
  cell_indx <- which(!is.na(cells))
  covsub <- cov[cell_indx, ]
  altsub <- alt[cell_indx, ]
  refsub <- ref[cell_indx, ]
  
  #summarising the matrices
  scov <- summary(covsub)
  salt <- summary(altsub)
  sref <- summary(refsub)
  
  #merging with i bring the cell and j being the position and x and y being the alt number and cov being total number (coverage)
  het_data <- merge(salt, scov, by = c("i", "j"), all = T) %>%
    # change NA values to 0
    replace_na(list(x.x = 0, x.y = 0)) %>%
    dplyr::rename(
      cell = i,
      snp_indx = j,
      altn = x.x,
      covn = x.y
    ) %>%
    merge(snps[, c(1, 2, 5)], by = 'snp_indx') %>%
    mutate(homozygosity = altn / covn) %>%
    merge(cn_df, by = 'cell') %>% #add column which is cov -alt
    mutate(refn = covn - altn)
  het_data$sample <- samp

  
  print("mut data gen")
  #somatic cells 
  somatic_cells <-metadata$cells[metadata$celltype %in% c("Early cyst", "Late cyst","Muscle") & metadata$sample == samp]
  #germ cells 
  germ_cells <- metadata$cells[!metadata$celltype %in% c("Early cyst", "Late cyst","Muscle") & metadata$sample == samp]
  
  
  #somatic genotype data 
  som_het_data <- filter(het_data, cell_barcode %in% somatic_cells) %>% 
    group_by(Pos, Chr) %>% 
    summarise(refn = sum(refn), 
              altn = sum(altn), 
              depth = sum(refn, altn)) %>% 
    filter(depth > 10) %>% 
    mutate(genotype = ifelse(refn == depth, "hom_ref", ifelse(altn == depth, "hom_alt", "het"))) %>% 
    filter(genotype != "het")
  
  
  #germ genotype data 
  germ_het_data <- filter(het_data, cell_barcode %in% germ_cells) %>% 
    merge(som_het_data[,c("Chr", "Pos", "genotype")], by = c("Chr", "Pos")) %>% 
    filter(covn > 1) %>% 
    filter(genotype != "het") %>% 
    mutate(cell_genotype = ifelse(refn == covn, "hom_ref", ifelse(altn == covn, "hom_alt", "het")))
  
  celltype_germ_data <- 
    germ_het_data %>% 
    mutate(Chr = ifelse(Chr %in% c("Chr_1", "Chr_2"), "A", "X")) %>% 
    group_by(cell_barcode, Chr) %>% 
    summarise(n_muts = sum(genotype != cell_genotype), 
              n_pos = n(), 
              mut_level = n_muts/n_pos) %>% 
    merge(metadata, by.x = 'cell_barcode', by.y = 'cells') 
  
  mut_data[[samp]] <- celltype_germ_data
  
  
}  


all_mut_data <- mut_data %>% 
  bind_rows()


all_mut_data %>% 
  mutate(celltype = factor(celltype, levels = c("GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  ggplot(aes(x = celltype, y = log(mut_level), fill = sample)) + 
  geom_boxplot()


all_mut_data %>% 
  mutate(celltype = factor(celltype, levels = c("GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  glmer(n_muts ~ celltype * Chr * treatment + (1|sample), poisson(link = "log"), data = ., 
      offset = log(n_pos)) %>% 
  summary()


load("data/trajectory/sce_GAMed_8k.RData")
traj_line <- SlingshotDataSet(sce)@curves$Lineage1
pts <- data.frame(Pseudotime = c(unlist(traj_line[[3]])), 
                  cells = names(c(unlist(traj_line[[3]]))))


Pseudotime_data <- reducedDims(sce)$UMAP %>% as.data.frame() %>% 
  merge(seurat_final@meta.data, by.x = 0, by.y = 'cells', all.y = F) %>% 
  mutate(cells = Row.names)


Pseudotime_data <- Pseudotime_data %>% 
  merge(pts, by = 'cells') %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate('Cell type' = celltype)


pt_amd <- all_mut_data %>% 
  merge(Pseudotime_data[,c("Pseudotime", "cells")], by.x = 'cell_barcode', by.y = 'cells')


pt_amd %>% 
  ggplot(aes(x = Pseudotime, y = log(mut_level))) + 
  geom_point()


pt_amd %>% 
  glmer(n_muts ~ Pseudotime * Chr + (1|treatment) + (1|sample), family = 'poisson', offset = log(n_pos), data = .) %>% 
  summary()


