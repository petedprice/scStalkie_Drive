
#cd /mnt/parscratch/users/bi1pp/Drive_pseudoreplicate
#singularity run ../Drive_var_rates/sing/rocker-tidyverse.img R
#install.packages("BiocManager", lib = "libs/")
library(BiocManager, lib = 'libs/')
#BiocManager::install("edgeR", lib = 'libs/')
library(edgeR, lib = 'libs/')
library(tidyverse)
rm(list = ls())

n=1000
replicates <- list()
counts <- read.table("pseudoreplicate/counts_pr.csv", sep = ',')
colnames(counts) <- gsub("\\.", "-", colnames(counts))
metadata <- read.table("pseudoreplicate/metadata_pr.csv", sep = ',')
metadata$celltype <- gsub("/", "and", metadata$celltype)
metadata$celltype <- gsub(" ", "_", metadata$celltype)

n=1000
i=1
ct='Muscle'
for (i in 1:n){
  print(i)
  for (ct in unique(metadata$celltype)){
    
    print(ct)
    cl_c=0
    cells_all <- list()
    for (trt in c("ST", "SR")){
      cl_c=cl_c+1
      cells <- c()
      for (s in 1:4){
        ps <- sample(metadata$cells[metadata$celltype == ct & metadata$treatment == "ST"], 100, replace = T)
        cells <-unique(c(cells,ps)) 
        
        s1 <- counts[,ps] %>% rowSums()
        if (s == 1 & trt == "ST"){
          s1_counts <- s1
        } else {
          s1_counts <- cbind(s1_counts, s1)
        }
        
      }
      cells_all[[cl_c]] <- cells
      
      
    }
    g1_keep <- (rowSums(counts[,cells_all[[1]]] > 1) > (0.05 * length(cells_all[[1]]))) #| 
    #rowSums(seurat_final@assays$RNA$counts[,ST_cells] > 1) > 10)
    
    g2_keep <- (rowSums(counts[,cells_all[[2]]] > 1) > (0.05 * length(cells_all[[2]])))# | 
    #                  rowSums(seurat_final@assays$RNA$counts[,SR_cells] > 1) > 10)
    
    keep <- g1_keep | g2_keep
    colnames(s1_counts) <- c(paste0("STa", 1:4), paste0("STb", 1:4))
    s1_counts <- s1_counts[keep,]
    y <- DGEList(s1_counts, samples = colnames(s1_counts), group = c(rep("STa", 4), rep("STb", 4)))
    y <- calcNormFactors(y)
    design <- model.matrix(~group, y$samples)
    y <- normLibSizes(y)
    y <- estimateDisp(y, design)
    cpm <- edgeR::cpm(y, log = F)
    
    keep1 <- cpm[,y$samples$group== "STa"] >= 1
    keep1 <- rowSums(keep1) > (ncol(keep1)/2)
    keep2 <- cpm[,y$samples$group== "STb"] >= 1
    keep2 <- rowSums(keep2) > (ncol(keep2)/2)
    keep <- keep1 | keep2
    
    y <- y[keep,keep.lib.sizes=FALSE]
    
    fit <- glmQLFit(y, design, robust = T)
    res <- glmQLFTest(fit, coef=ncol(design))
    toptags <- topTags(res, n = nrow(res))
    out <- toptags$table %>% 
      mutate(gene = rownames(.), 
             celltype = ct, 
             rep = paste0("rep", i), 
             treatment = "STaSTb", 
             time = date())
    
    filename = paste0("outdata/", gsub(" ", "", date()), 
                      "rep", i, ct, ".csv")
    
    write.table(out, filename, sep = ',', quote = F)
  }

}




#### ANALYSING LOCALLY 
library(BiocManager, lib = 'libs/')
library(edgeR, lib = 'libs/')
library(tidyverse)
rm(list = ls())

data_files <- list.files("outdata/", full.names = T, pattern = ".csv")

ortholog_table <- read.table("pseudoreplicate/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "") %>% 
  dplyr::select(REF_GENE_NAME, chr) %>% unique() %>% 
  mutate(REF_GENE_NAME = gsub("_", "-", REF_GENE_NAME))

sum_data_func <- function(file, ot = ortholog_table){
  tmp <- read.csv(file) %>% 
    merge(ot, by.x = 'gene', by.y = 'REF_GENE_NAME') %>% 
    mutate(sig = ifelse(logFC >= 1 & FDR < 0.05, "siga",
                        ifelse(logFC <=-1 & FDR < 0.05, "sigb", "non"))) %>%
    group_by(celltype, rep, sig, chr, treatment, time) %>%
    summarise(n = n()) %>%
    group_by(celltype, rep, treatment, chr, time) %>%
    summarise(prop = n/sum(n),
              sig = sig,
              n = n)
  filename = paste0("sumarised/", tmp$celltype[1], tmp$rep[1], tmp$treatment[1], gsub(" ", "", date()), "_sum_data.csv")
  write.csv(tmp, filename, quote = F, row.names = F)
  system(paste0("mv ", file, " outdata/completed"))
}

lapply(data_files, sum_data_func, ot = ortholog_table) 


summarys <- lapply(list.files("data/pseudoreplicate/sumarised/", full.names = T), read.csv) %>% 
  bind_rows()


comparisons <- list(c("STvsSR", "STvsST"), c("Autosomes", "X_chromosome"))
summarys %>% dim()
table(summarys$treatment, summarys$celltype)


summarys %>% 
  group_by(celltype, treatment, chr) %>% 
  summarise(maxrep = max(rep),
            reps = length((rep))) %>% View()


comparisons1 <- list(c("STvsST", "STvsSR"))
comparisons2 <- list(c("X_chromosome", "Autosomes"))


  
summarys <- summarys %>% 
  mutate(treatment = ifelse(treatment == "STaSTb", "STvsST", "STvsSR"),
         celltype = ifelse(celltype == "GSCandSpermatogonia", "GSC/Spermatogonia", celltype),
         celltype = gsub("_", " ", celltype),
         celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids")),
         chr = ifelse(chr != "Chr_X", "Autosomes", "X_chromosome"))   


summarys %>% filter(sig != "non") %>% 
  ggplot(aes(x = treatment, y = prop * 100, fill = chr)) + 
  geom_boxplot(outlier.alpha = 0.011, position = position_dodge(0.75)) +
  scale_fill_grey(start = 0.4, end = 0.7) + 
  theme_classic() + 
  facet_wrap(~celltype) + 
  labs(fill = "", y = "% of genes differentially expressed") + 
  
  # Wilcoxon between treatment groups (STvsST vs STvsSR), within each chromosome
  stat_compare_means(method = "wilcox.test", 
                    comparisons = comparisons1, 
                    aes(label = after_stat(p.signif), group = treatment), 
                    label.y = 32,
                    position = position_dodge(0.4)) +
  stat_compare_means(method = "wilcox.test", 
                    aes(label = after_stat(p.signif)), 
                    label.y = 28) + 
  ylim(0,45)
                  



ggsave("plots/manuscript_plots/SX_pseudoreplicate.pdf", height = 7, width = 7)
system("open plots/manuscript_plots/SX_pseudoreplicate.pdf")




