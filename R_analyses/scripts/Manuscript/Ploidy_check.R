library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(gridExtra)
library(cowplot)
rm(list = ls())
load("data/RData/seurat_final.RData")
samples <- c("sr1", "sr2", "sr3", "sr5", "st1","st2","st3","st5")

run_checks = "yes"
if (run_checks == "yes"){
  
  metadata <- seurat_final@meta.data
  thresholds <- data.frame(het = c(0,0,0,1,1), hom = c(0,1,1,3,3), nsnps = rep(c(0,5,10,5,10), each = 1))
  metadatas <- list()
  het_data_sum_save <- list()
  het_data_save <- list()
  c=0
  #samp = samples[1]
  for (samp in samples){  
    print(samp)
    snps <- read.table(paste0("data/scAlleleCount/relaxed/data/", samp, "_filt.recode.bed"))
    snps$snp_indx <- 1:nrow(snps)
    colnames(snps) <- c("Chr", "Pos", "Ref", "Alt", "snp_indx")
    barcodes <- read.table(paste0("data/scAlleleCount/relaxed///data/", samp, "_barcodes.tsv"))[,1]
    cov <- readMM(file = paste0("data/scAlleleCount/relaxed///data/", samp, "covmat.mtx"))
    alt <- readMM(file = paste0("data/scAlleleCount/relaxed//data/", samp, "altmat.mtx"))
    ref <- readMM(file = paste0("data/scAlleleCount/relaxed//data/", samp, "refmat.mtx"))
    
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
          hetn = length(which(homozygosity > 0 & homozygosity < 1 & altn > thresholds[tt,1] & refn > thresholds[tt,1])), #& covn > thresholds[tt,2])),
          hetd= sum(covn[!homozygosity %in% c(0,1) & covn > thresholds[tt,2]]),
          total_hethoms = (hetn + homn), 
          hethomd=homd+hetd,
          readsn = sum(covn),
          snpsn = length(unique(snp_indx)),
          snpcov = hethomd/total_hethoms) %>% 
        mutate(heterozygosity = homn/total_hethoms)
      het_data_sum$coverage_check <- "OK"
      het_data_sum$coverage_check[het_data_sum$total_hethoms < thresholds[tt,3] | is.na(het_data_sum$snpsn)] <- "low_coverage"
      het_data_sum$sample <- samp
      het_data_sum$het_threshold <- thresholds[tt,1] + 1
      het_data_sum$hom_threshold <- thresholds[tt,2] + 1
      het_data_sum$snp_threshold <- thresholds[tt,3]
      
      
      het_data_sum_save[[c]] <- het_data_sum
    }
  }
    
    
    
  new_meta_data <- bind_rows(het_data_sum_save)
  new_meta_data <- new_meta_data %>%  
    merge(metadata, by.x = c('cell_barcode', 'sample'), by.y = c('cells', 'sample'), all.y = F)
  new_meta_data$ploidy_class <- ifelse(new_meta_data$heterozygosity < 0.95, "diploid", "haploid")
  new_meta_data$ploidy_class[is.na(new_meta_data$ploidy_class)] <- "Missing"
  
  summary_table <- new_meta_data %>% 
    group_by(het_threshold, hom_threshold, snp_threshold, celltype, ploidy_class) %>%
    summarise(n =n())
  grid.table(summary_table)
  
  save(summary_table, new_meta_data, file = "data/scAlleleCount//ploidy_params_check.RData")
} else if (run_checks == "no"){
  load(file = "data/scAlleleCount//ploidy_params_check.RData")
}



new_meta_data$celltype <- factor(new_meta_data$celltype, levels = 
                                   c("Muscle", "Early cyst", "Late cyst", "GSC/Spermatogonia", 
                                     "Primary spermatocytes", "Secondary spermatocytes", "Early spermatids", "Late spermatids"))
###########################

p1 <- 
  new_meta_data %>% 
  filter(treatment == "ST") %>% 
  filter(snp_threshold ==0 & het_threshold ==1 & hom_threshold == 1) %>% 
  filter(Chr != "Chr_X") %>% 
  mutate(Ploidy = ploidy_class) %>% 
  filter(celltype %in% c("Early cyst", "Late cyst", "Muscle", "GSC/Spermatogonia", "Primary spermatocytes")) %>% 
  ggplot(aes(x = snpcov, y = 100*heterozygosity, colour = Ploidy, fill = Ploidy)) +
  geom_point(alpha = 0.5) + 
  scale_colour_brewer(palette = 'Set1') + #, guide = 'none') +
  scale_fill_brewer(palette = 'Set1',guide = guide_legend(override.aes = list(linetype = c(0, 0),
                                                                              shape = c(NA, NA),
                                                                              color = c(NA, NA), 
                                                                              alpha = c(1,1)))) +
  #geom_smooth(method = "lm", colour = 'red', aes(group = 1)) + 
  stat_cor(method = 'pearson', aes(group = 1), label.y = 87) + 
  labs(y = "% of homozygous sites per cell", x = "Mean per-cell SNP coverage") + 
  theme_classic() +  
  theme(legend.position = c(0.8,0.8))

p2 <- new_meta_data %>% 
  filter(treatment == "ST") %>% 
  filter(snp_threshold ==0 & het_threshold ==1 & hom_threshold == 1) %>% 
  mutate(Ploidy = ploidy_class) %>% 
  filter(Chr != "Chr_X") %>% 
  #filter(celltype %in% c("Early cyst", "Late cyst", "Muscle", "GSC/Spermatogonia", "Primary spermatocytes")) %>% 
  ggplot(aes(x = celltype, fill = Ploidy)) + 
  geom_bar(position = 'dodge') + 
  scale_fill_brewer(palette = 'Set1') + 
  theme_classic() + 
  labs(y = "Number of cells", x = "") + 
  theme(legend.position = 'none', 
        axis.text.x=element_blank())

p3 <- new_meta_data %>% 
  filter(treatment == "ST") %>% 
  filter(ploidy_class != "Missing") %>%
  filter(snp_threshold ==10 & het_threshold == 2& hom_threshold == 4) %>% 
  filter(Chr != "Chr_X") %>% 
  mutate(Ploidy = ploidy_class) %>% 
  filter(celltype %in% c("Early cyst", "Late cyst", "Muscle", "GSC/Spermatogonia", "Primary spermatocytes")) %>% 
  ggplot(aes(x = snpcov, y = 100*heterozygosity, colour = Ploidy)) +
  geom_point(alpha = 0.2) + 
  scale_colour_brewer(palette = 'Set1') +
  #geom_smooth(method = "loess", colour = 'red', aes(group = 1)) + 
  stat_cor(method = 'pearson', aes(group = 1), label.y = 5) + 
  labs(y = "% of homozygous sites  per cell", x = "Mean per-cell SNP coverage") + 
  theme_classic() + 
  theme(legend.position = 'none')

p4 <- new_meta_data %>% 
  filter(treatment == "ST") %>% 
  filter(ploidy_class != "Missing") %>%
  filter(snp_threshold ==10 & het_threshold == 2& hom_threshold == 4) %>% 
  mutate(Ploidy = ploidy_class) %>% 
  filter(Chr != "Chr_X") %>% 
  #filter(celltype %in% c("Early cyst", "Late cyst", "Muscle", "GSC/Spermatogonia", "Primary spermatocytes")) %>% 
  ggplot(aes(x = celltype, fill = Ploidy)) + 
  geom_bar(position = 'dodge') + 
  scale_fill_brewer(palette = 'Set1') + 
  theme_classic() + 
  labs(y = "Number of cells", x = "") + 
  guides(fill = guide_legend(title = "Predicted ploidy")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  theme(legend.position = 'none')

new_meta_data %>% 
  filter(treatment == "ST") %>% 
  filter(ploidy_class != "Missing") %>%
  filter(snp_threshold ==10 & het_threshold == 2& hom_threshold == 4) %>% 
  mutate(Ploidy = ploidy_class) %>% 
  filter(Chr != "Chr_X") %>%  dim()
new_meta_data %>% 
  filter(treatment == "ST") %>% 
  filter(snp_threshold ==0 & het_threshold ==1 & hom_threshold == 1) %>% 
  mutate(Ploidy = ploidy_class) %>% 
  filter(Chr != "Chr_X") %>% dim()

ab <- align_plots(p1,p2, align = 'h', axis = 'b')
cd <- align_plots(p3,p4, align = 'h', axis = 'b')
bd <- plot_grid(ab[[2]],cd[[2]], labels = c("(b)", "(d)"), align = 'v', 
                axis = 'r', ncol = 1, label_x = -0.03, rel_heights = c(4,5))
ac <- plot_grid(ab[[1]],cd[[1]], labels = c("(a)", "(c)"), align = 'v', 
                axis = 'r', ncol = 1, label_x = -0.03, rel_heights = c(4,5))
Ploidy_Sup <- plot_grid(ac, bd, rel_widths = c(4,4))

ggsave("plots/Ploidy_Sup.svg", Ploidy_Sup, width = 8, height = 9)
system("open plots/Ploidy_Sup.svg")

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

new_meta_data %>% 
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



