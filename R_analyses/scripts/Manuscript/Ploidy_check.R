library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(gridExtra)
rm(list = ls())
load("data/RData/seurat_final.RData")
samples <- c("sr1", "sr2", "sr3", "sr5", "st1","st2","st3","st5")

run_checks = "yes"
if (run_checks == "yes"){
  
  metadata <- seurat_final@meta.data
  #metadata$ploidy <- NA
  #metadata$nsnps <- NA
  
  #thresholds <- data.frame(het = c(0,1,2,3), hom = c(0,2,4,6), nsnps = rep(c(0,5,10), each = 4))
  thresholds <- data.frame(het = c(0), hom = c(0), nsnps = rep(c(0), each = 1))
  metadatas <- list()
  het_data_sum_save <- list()
  het_data_save <- list()
  c=0
  #samp = samples[1]
  for (samp in samples){  
    print(samp)
    snps <- read.table(paste0("data/scAlleleCount/relaxed/", samp, "_snps.txt"))
    snps$snp_indx <- 1:nrow(snps)
    colnames(snps) <- c("Chr", "Pos", "Ref", "Alt", "snp_indx")
    barcodes <- read.table(paste0("data/scAlleleCount/relaxed//data/", samp, "_barcodes.tsv"))[,1]
    cov <- readMM(file = paste0("data/scAlleleCount/relaxed/data/", samp, "covmat.mtx"))
    alt <- readMM(file = paste0("data/scAlleleCount/relaxed/data/", samp, "altmat.mtx"))
    ref <- readMM(file = paste0("data/scAlleleCount/relaxed/data/", samp, "refmat.mtx"))
    
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
    scov <- summary(covsub)
    salt <- summary(altsub)
    sref <- summary(refsub)
    
    #merging with i bring the cell and j being the position and x and y being the alt number and cov being total number (coverage)
    het_data <- merge(salt, scov, by=c("i", "j"), all = T) %>% 
      # change NA values to 0
      replace_na(list(x.x = 0, x.y = 0)) %>% 
      dplyr::rename(cell = i, snp_indx = j, altn = x.x, covn = x.y) %>%
      merge(snps[,c(1,2,5)], by = 'snp_indx') %>% 
      mutate(homozygosity = altn/covn) %>% 
      merge(cn_df, by = 'cell') %>% #add column which is cov -alt 
      mutate(refn = covn - altn)
    het_data$sample <- samp
    het_data_save[[samp]] <- het_data
    for (tt in 1:nrow(thresholds)){
      print(tt)
      c=c+1
      het_data_sum <- het_data %>% 
        mutate(Chr = ifelse(Chr == "Chr_X", "Chr_X", "Chr_A")) %>% 
        group_by(cell_barcode, Chr) %>% 
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
    
    
    
  new_meta_data <- bind_rows(het_data_sum_save)
  new_meta_data <- new_meta_data %>%  
    merge(metadata, by.x = c('cell_barcode', 'sample'), by.y = c('cells', 'sample'), all.y = F)
  new_meta_data$ploidy_class <- ifelse(new_meta_data$heterozygosity < 0.99, "diploid", "haploid")
  new_meta_data$ploidy_class[is.na(new_meta_data$ploidy_class)] <- "Missing"
  
  summary_table <- new_meta_data %>% 
    group_by(depth_threshold, snp_threshold, celltype, ploidy_class) %>%
    summarise(n =n())
  grid.table(summary_table)
  
  save(summary_table, new_meta_data, file = "data/scAlleleCount//ploidy_params_check.RData")
} else if (run_checks == "no"){
  load(file = "data/scAlleleCount//ploidy_params_check.RData")
}

bar_graph <- summary_table %>% 
  ggplot(aes(x = celltype, fill = ploidy_class, y = n)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~depth_threshold + snp_threshold, nrow = 4)
st_ss <- filter(summary_table, snp_threshold %in% c(5,10) & depth_threshold %in% c(c,1,2, 3)) %>% 
  spread(ploidy_class, n) %>% 
  mutate(total = sum(diploid, haploid, Missing), 
         prop_diploid = diploid/total,
         prop_haploid = haploid/total, 
         prop_missing = Missing/total)
row.names(st_ss) <- NULL
png("plots/scAlleleCount_params_table.png", height = 30*nrow(st_ss), width = 120*ncol(st_ss))
grid.table(st_ss, rows = NULL)
dev.off()


  
ggsave("plots/ploidy_check.png", bar_graph, width = 40, height = 20, units = "in")

write.table(summary_table, 'data/scAlleleCount/ploidy_check_summary.tsv', sep = '\t', quote = F, row.names = F)


###########################

new_meta_data %>% 
  filter(Chr != "Chr_X") %>% 
  filter(celltype %in% c("Early cyst", "Late cyst", "Muscle")) %>% 
  ggplot(aes(x = snpcov, y = heterozygosity)) +
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm", colour = 'red') + 
  stat_cor(method = 'spearman') + 
  labs(y = "Heterozygosity", x = "Mean per-cell SNP coverage")



############################
new_meta_data %>% 
  filter(Chr == "Chr_X") %>% 
  filter(depth_threshold %in% c(0) & snp_threshold %in% c(0)) %>% 
  ggplot(aes(x = celltype, y = heterozygosity, fill = treatment)) + 
  geom_boxplot() + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 1.01, method = "wilcox.test")
                      
             
new_meta_data %>% 
  filter(Chr == "Chr_X") %>% 
  filter(depth_threshold %in% c(0) & snp_threshold %in% c(0)) %>% 
  ggplot(aes(x = -log(heterozygosity))) + geom_density()
############################

################## JO COPY GRAPH ########################

new_meta_data %>% filter(depth_threshold %in% c(0,3) & snp_threshold %in% c(0)) %>% 
  filter(coverage_check == "OK") %>% 
  ggplot(aes(x = log(hethomd), y = heterozygosity)) + 
  geom_point(alpha = 0.4) + 
  geom_smooth(method = 'lm') + #stat_cor(method = "spearman", position = 'jitter') + 
  ggpmisc::stat_correlation(method = 'pearson', aes(
    label = paste(after_stat(r.label),
                  after_stat(p.value.label),
                  after_stat(n.label),
                  sep = "*\", \"*")), position = 'jitter') +
  labs(y = "homozygosity, 1=all sites hom, 0=all sites het",
       x = "log(number of snps/summed depth over those snps)") + 
  facet_wrap(~depth_threshold + snp_threshold, scales = "free_y")
  

#################################



