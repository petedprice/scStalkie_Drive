library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
load('data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned_celltype.RData')

samples <- c("sr1", "sr2", "sr3", "sr5", "st1","st2","st3","st5")
metadata <- seurat_integrated@meta.data
metadata$ploidy <- NA
metadata$nsnps <- NA
het_data_sum_save <- list()
het_data_save <- list()

for (samp in samples){  
  print(samp)
  snps <- read.table(paste0("data/scAlleleCount/stringent/", samp, "_snps.txt"))
  snps$snp_indx <- 1:nrow(snps)
  colnames(snps) <- c("Chr", "Pos", "Ref", "Alt", "snp_indx")
  barcodes <- read.table(paste0("data/scAlleleCount/stringent//data/", samp, "_barcodes.tsv"))[,1]
  cov <- readMM(file = paste0("data/scAlleleCount/stringent/data/", samp, "covmat.mtx"))
  alt <- readMM(file = paste0("data/scAlleleCount/stringent/data/", samp, "altmat.mtx"))
  ref <- readMM(file = paste0("data/scAlleleCount/stringent/data/", samp, "refmat.mtx"))
  
  #subsetting matrices for cells we have present in our seurat object
  cells <- match(paste0(samp, "_", barcodes), rownames(seurat_integrated@meta.data))
  cell_names <- rownames(seurat_integrated@meta.data)[cells]
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
    rename(cell = i, snp_indx = j, alt = x.x, cov = x.y) %>%
    merge(snps[,c(1,2,5)], by = 'snp_indx') %>% 
    mutate(homozygosity = alt/cov) %>% 
    filter(cov > 3) %>% 
    merge(cn_df, by = 'cell') %>% #add column which is cov -alt 
    mutate(ref = cov - alt)
  het_data$sample <- samp
  het_data_save[[samp]] <- het_data
  
  het_data_sum <- het_data %>% 
    group_by(cell_barcode) %>% 
    summarise( #summarise so hom is no of pos where homozygosity is 1 and het is no of pos where homozygosity is not 1
      hom = length(which(homozygosity %in% c(0,1))), 
      het = length(which(homozygosity > 0 & homozygosity < 1 & alt > 1 & ref > 1)),
      total_hethoms = (het + hom), 
      snp_cov = sum(cov)/length(snp_indx), 
      total_cov = sum(cov),
      nsnps = length(snp_indx)) %>% 
    mutate(ploidy = hom/total_hethoms)
  het_data_sum$sample <- samp
  het_data_sum$ploidy_class <- ifelse(het_data_sum$ploidy < 0.95, "diploid", "haploid")
  het_data_sum_save[[samp]] <- het_data_sum
}



new_meta_data <- bind_rows(het_data_sum_save) %>% 
  merge(seurat_integrated@meta.data, by.x = 'cell_barcode', by.y = 'cells', all.y = T)

new_meta_data$coverage <- "OK"
new_meta_data$coverage[new_meta_data$nsnps < 10 | is.na(new_meta_data$nsnps)] <- "low_coverage"
new_meta_data$ploidy_class <- ifelse(new_meta_data$ploidy < 0.95, "diploid", "haploid")
new_meta_data$ploidy_class[is.na(new_meta_data$ploidy_class)] <- "Missing"

rownames(new_meta_data) <- new_meta_data$cell_barcode
si_ploidy <- seurat_integrated 
si_ploidy <- AddMetaData(si_ploidy, new_meta_data)
plots <- list()
plots[[1]] <- DimPlot(si_ploidy, group.by = 'ploidy_class')
plots[[2]] <- si_ploidy@meta.data %>% 
  ggplot(aes(x = ploidy_class, fill = coverage)) + geom_bar(position = 'dodge') + theme_minimal() +
  facet_wrap(~celltype)
plots[[3]] <- si_ploidy@meta.data %>% 
  filter(!is.na(snp_cov)) %>% 
  ggplot(aes(x = (snp_cov), colour = ploidy_class)) + geom_histogram() +
  facet_wrap(~ploidy_class)

plots %>% 
  ggarrange(plotlist = .) %>% 
  ggsave("plots/ploidy_check_stringent.pdf", .)

seurat_integrated@meta.data <- new_meta_data
DimPlot(seurat_integrated, group.by = 'ploidy_class')
bin_data <- list()
for (ct in unique(new_meta_data$celltype)){
  bin_df <- data.frame(depth = seq(0,1000,5))
  bin_df$dip_prob <- sapply(bin_df$depth, 
                      function(x)
                        nrow(filter(new_meta_data, 
                                    total_cov > x & total_cov < x+5 & ploidy_class == 'diploid'))/nrow(filter(new_meta_data, total_cov > x & total_cov < x+5)))
  bin_df$ncells <- sapply(bin_df$depth,
                          function(x)
                            nrow(filter(new_meta_data, 
                                        total_cov > x & total_cov < x+5 & celltype == ct)))
  bin_df$celltype <- ct
  bin_data[[ct]] <- bin_df
}

bin_data_df <- bin_data %>% bind_rows()
plots <- list()
plots[[1]] <- bin_data_df %>% 
  filter(ncells > 1) %>% 
  ggplot(aes(x = depth, y = dip_prob, colour = log(ncells)))+ geom_point() + 
  facet_wrap(~celltype) + #add new colour scale for ncells that goes across colours 
  scale_colour_gradient(low = "blue", high = "orange")
plots[[2]] <- new_meta_data %>% 
  ggplot(aes(x = celltype, y = log(nsnps))) + 
  geom_boxplot()
plots %>% 
  ggarrange(plotlist = ., nrow = 2) %>% 
  ggsave("plots/coverage_ploidy.pdf", .)


#Distribution of coverage for the SNPs used to call diploid/haploid in each cell type?
all_snp_data <- bind_rows(het_data_save)
ploidy_vs_cov_plot <- list()
ploidy_vs_cov_plot[[1]] <- all_snp_data %>% 
  reframe(zygosity = homozygosity, 
          coverage = cov, zygosity = ifelse(zygosity > 0.5, 1-zygosity, zygosity)) %>% 
  ggplot(aes(x = as.factor(coverage), y = zygosity)) + geom_boxplot() + 
  labs(x = 'per site read coverage', y = 'zygosity')



ploidy_vs_cov_plot[[2]] <- new_meta_data %>% 
  ggplot(aes(x = log(total_hethoms), y = ploidy, colour = ploidy_class)) + geom_point(alpha = 0.1) + 
  labs(x = "log(number of positions)", y = "ploidy") + facet_wrap(~celltype)

ploidy_vs_cov_plot[[3]] <- new_meta_data %>% 
  ggplot(aes(x = log(total_hethoms), fill = ploidy_class)) + geom_density(alpha = 0.5) +
  labs(x = "log(number of positions)") +
  facet_wrap(~celltype)

ploidy_vs_cov_plot %>%
  ggarrange(plotlist = ., nrow = 10) %>% 
  ggsave("plots/ploidy_vs_cov.png", ., width = 20, height = 40)
