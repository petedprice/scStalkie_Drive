
#cd /mnt/parscratch/users/bi1pp/Drive_pseudoreplicate
#singularity run ../Drive_var_rates/sing/rocker-tidyverse.img R
#install.packages("BiocManager", lib = "libs/")
library(BiocManager, lib = 'libs/')
#BiocManager::install("edgeR", lib = 'libs/')
library(edgeR, lib = 'libs/')
library(tidyverse)
rm(list = ls())

counts <- read.table("pseudoreplicate/counts_pr.csv", sep = ',')
colnames(counts) <- gsub("\\.", "-", colnames(counts))
metadata <- read.table("pseudoreplicate/metadata_pr.csv", sep = ',')
metadata$celltype <- gsub("/", "and", metadata$celltype)
metadata$celltype <- gsub(" ", "_", metadata$celltype)
ortholog_table <- read.table("pseudoreplicate/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "") %>% 
  mutate(chr = ifelse(chr == "Chr_X", "Chr_X", "Chr_A")) %>% 
  dplyr::select(REF_GENE_NAME, chr) %>% unique() %>% 
  mutate(REF_GENE_NAME = gsub("_", "-", REF_GENE_NAME))

n=1000
#for (comp in list(c("ST", "SR"), c("ST", "ST"))){
for (comp in list(c("ST", "ST"))){
  for (i in 1:n){
    print(comp)
    print(i)
    for (ct in unique(metadata$celltype)){
      
      print(ct)
      cl_c=0
      cells_all <- list()
      for (trt in 1:length(comp)){
        cl_c=cl_c+1
        cells <- c()
        for (s in 1:4){
          ps <- sample(metadata$cells[metadata$celltype == ct & metadata$treatment == comp[trt]], 100, replace = T)
          cells <-unique(c(cells,ps)) 
          
          s1 <- counts[,ps] %>% rowSums()
          if (cl_c == 1 & s ==1){
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
      names <- c(paste0(comp[1],"a", 1:4), paste0(comp[2],"b", 1:4))
      groups <- c(rep(paste0(comp[1], "a"), 4), rep(paste0(comp[2], "b"), 4))
      colnames(s1_counts) <- names
      s1_counts <- s1_counts[keep,]
      y <- DGEList(s1_counts, samples = colnames(s1_counts), group = groups)
      y <- calcNormFactors(y)
      design <- model.matrix(~group, y$samples)
      y <- normLibSizes(y)
      y <- estimateDisp(y, design)
      cpm <- edgeR::cpm(y, log = F)
      
      keep1 <- cpm[,y$samples$group== paste0(comp[1], "a")] >= 1
      keep1 <- rowSums(keep1) > (ncol(keep1)/2)
      keep2 <- cpm[,y$samples$group== paste0(comp[2], "b")] >= 1
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
               treatment = paste(comp, collapse = "vs")) %>% 
        merge(ortholog_table, by.x = 'gene', by.y = 'REF_GENE_NAME') %>% 
        mutate(sig = ifelse(logFC >= 1 & FDR < 0.05, "siga",
                            ifelse(logFC <=-1 & FDR < 0.05, "sigb", "non"))) %>%
        group_by(celltype, rep, sig, chr, treatment) %>%
        reframe(n = n()) %>%
        group_by(celltype, rep, treatment, chr) %>%
        reframe(prop = n/sum(n),
                  sig = sig,
                  n = n)
        

      filename = paste0("sumarised/", out$celltype[1], out$rep[1], out$treatment[1], gsub(" ", "", date()), "_sum_data.csv")
      write.csv(out, filename, quote = F, row.names = F)
    }
  }
}



