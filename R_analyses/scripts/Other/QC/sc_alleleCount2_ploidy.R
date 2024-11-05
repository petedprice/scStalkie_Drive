library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(gridExtra)

load('data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned_celltype.RData')

samples <- c("sr1", "sr2", "sr3", "sr5", "st1","st2","st3","st5")

run_checks = "no"
if (run_checks == "yes"){
  metadata <- seurat_integrated@meta.data
  #metadata$ploidy <- NA
  #metadata$nsnps <- NA
  
  thresholds <- data.frame(het = c(0,1,2,3), hom = c(0,2,4,6), nsnps = rep(c(0,5,10), each = 4))
  metadatas <- list()
  het_data_sum_save <- list()
  het_data_save <- list()
  c=0
  #samp = samples[1]
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
      rename(cell = i, snp_indx = j, altn = x.x, covn = x.y) %>%
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
        group_by(cell_barcode) %>% 
        summarise( #summarise so hom is no of pos where homozygosity is 1 and het is no of pos where homozygosity is not 1
          homn = length(which(homozygosity %in% c(0,1) & covn > thresholds[tt,2])), 
          homd= sum(covn[homozygosity %in% c(0,1) & covn > thresholds[tt,2]]),
          hetn = length(which(homozygosity > 0 & homozygosity < 1 & altn > thresholds[tt,1] & refn > thresholds[tt,1] & covn > thresholds[tt,2])),
          hetd= sum(covn[homozygosity %in% c(0,1) & covn > thresholds[tt,2]]),
          total_hethoms = (hetn + homn), 
          hethomd=homd+hetd,
          snpcov = sum(covn)/length(total_hethoms), 
          readsn = sum(covn),
          snpsn = length(unique(snp_indx))) %>% 
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
  new_meta_data$ploidy_class <- ifelse(new_meta_data$heterozygosity < 0.95, "diploid", "haploid")
  new_meta_data$ploidy_class[is.na(new_meta_data$ploidy_class)] <- "Missing"
  
  summary_table <- new_meta_data %>% 
    group_by(depth_threshold, snp_threshold, celltype, ploidy_class) %>%
    summarise(n =n())
  grid.table(summary_table)
  
  save(summary_table, new_meta_data, file = "data/scAlleleCount//ploidy_params_check.RData")
}



load(file = "data/scAlleleCount//ploidy_params_check.RData")

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



################## JO COPY GRAPH ########################

new_meta_data %>% filter(depth_threshold %in% c(0,3) & snp_threshold %in% c(5,10)) %>% 
  filter(coverage_check == "OK") %>% 
  ggplot(aes(x = log(hethomd), y = heterozygosity)) + 
  geom_point(alpha = 0.4) + 
  geom_smooth(method = 'lm') + #stat_cor(method = "spearman", position = 'jitter') + 
  stat_correlation(method = 'pearson', aes(
    label = paste(after_stat(r.label),
                  after_stat(p.value.label),
                  after_stat(n.label),
                  sep = "*\", \"*")), position = 'jitter') +
  labs(y = "homozygosity, 1=all sites hom, 0=all sites het",
       x = "log(number of snps/summed depth over those snps)") + 
  facet_wrap(~depth_threshold + snp_threshold, scales = "free_y")
  

#################################



new_meta_data %>% 
  filter(snp_threshold == 0, depth_threshold == 0) %>% 
  summarise(min = min(nsnps))
  group_by(snp_threshold) %>% 
  summarise(n = length(which(ploidy_class == 'diploid')))

################### OLD RANDOM PLOTS ----------------------

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
  geom_boxplot() + 
  facet_wrap(~coverage)
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




###########


all_snp_data_celltype <- all_snp_data %>% 
  merge(seurat_integrated@meta.data, by.x = 'cell_barcode', by.y = 'cell_barcode', all.y = T)
stringency <- list()
for (i in c(1,2,3,4,5)){
  all_snp_data_celltype$hethom <- NA
  all_snp_data_celltype$hethom[all_snp_data_celltype$homozygosity %in% c(0,1) & all_snp_data_celltype$cov > 1] <- "hom"
  all_snp_data_celltype$hethom[!all_snp_data_celltype$homozygosity %in% c(0,1) & all_snp_data_celltype$ref > 2 & all_snp_data_celltype$alt > 2] <- "het"



  all_snp_data_celltype %>% 
    group_by(celltype, hethom) %>% 
    summarise(n = n()) %>%  View()
}

asdc_plot <- all_snp_data_celltype %>% 
  #filter(!is.na(hethom)) %>% 
  ggplot(aes(x = log(cov), colour = hethom)) + geom_density() + 
  facet_wrap(~celltype)

asdc_plot <- all_snp_data_celltype %>% 
  ggplot(aes(x = log(cov), y = homozygosity)) + geom_point() + 
  facet_wrap(~celltype)


all_snp_data_celltype$hethom <- NA
all_snp_data_celltype$hethom[all_snp_data_celltype$homozygosity %in% c(0,1) & all_snp_data_celltype$cov > 3] <- "hom"
all_snp_data_celltype$hethom[!all_snp_data_celltype$homozygosity %in% c(0,1) & all_snp_data_celltype$ref > 3 & all_snp_data_celltype$alt > 3] <- "het"

all_snp_data_celltype %>% 
  group


all_snp_data_celltype %>% 
  group_by(celltype, hethom) %>% 
  summarise(n = n()) %>%  View()


################## number of snps per celltype distribution ------------
new_meta_data %>% 
  ggplot(aes(x = log(total_hethoms), colour = ploidy_class)) + geom_density() +
  facet_wrap(~celltype)


new_meta_data %>% 
  ggplot(aes(x = celltype, colour = ploidy_class)) + geom_bar()

new_meta_data %>% 
  ggplot(aes(x = celltype, y = ploidy)) + geom_boxplot() 



