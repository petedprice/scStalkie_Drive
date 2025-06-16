################# load libraries and functions and data prep ###################
library(dplyr)
library(Seurat)
library(openxlsx)
library(ggpubr)
library(scuttle)
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(ggrepel)
library(GenomicFeatures)
library(ggrepel)
library(cowplot)
rm(list = ls())

#### LOAD DATA AND PREP DATA ---
load("data/RData/seurat_final.RData")
load("data/RData/DEG_DC.RData")

seurat_final <- JoinLayers(seurat_final)
Idents(seurat_final) <- "celltype"
levels(seurat_final) <- c("Muscle", "Early cyst", "Late cyst", 
                          "GSC/Spermatogonia", "Primary spermatocytes", 
                          "Secondary spermatocytes", "Early spermatids", "Late spermatids")
seurat_final@meta.data <- seurat_final@meta.data %>%  
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids")))

DefaultAssay(seurat_final) <- "RNA"

ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")


ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]


genes_search <- c("g1937", "PB.151","gene_10203")
ortholog_table %>% 
  filter(REF_GENE_NAME %in% genes_search) %>% 
  group_by(REF_GENE_NAME) %>% 
  summarise(start = min(start), 
            end = max(end)) %>% 
  mutate(search_term = 
           paste0("Chr_X:", start, "-", end))

tidy_cpm_edger %>% 
  filter(genes %in% genes_search)


DGE_inversion <- dif_exp_data

new_inversions <- data.frame(inversions = c("Tdal_INV0", "Tdal_INV1", "Tdal_INV2", "Tdal_INV3", "Tdal_INV_TLOC"), 
                             start = c(0, 398747, 32269986, 59107670, 76695356), 
                             end = c(398747, 9036332, 48326019, 87076097, 80442481), 
                             start2 = c(17545, 1252862, 32269991, 59117787, 76695380), 
                             end2 = c(1252862, 11084918, 48326025, 87076515, 80443921)) %>% 
  mutate(start_uncert = start2-start, 
         end_uncert = end2-end)
a <- new_inversions[,c(1,2,4)]
b <- new_inversions[,c(1,3,5)]
colnames(a) <- c("Inversion", "start", "end")
colnames(b) <- c("Inversion", "start", "end")

new_invs2 <- rbind(a,b)
new_invs2$type <- c(rep("Start", 5), rep("Finish", 5))
dif_exp_data %>% 
  filter(chr == "Chr_X") %>% 
  apply(., 1, function(x){
    inv_points <- dplyr::filter(new_invs2, 
                         start < x[17] & end > x[17] | 
                           start < x[18] & end > x[18] | 
                           start < x[17] & end > x[18] | 
                           start > x[17] & end < x[18])
    if (nrow(inv_points) > 0){
      out <- cbind(x, inv_points)
      return(out)
    }
  }
  )
data <- dif_exp_data[dif_exp_data$chr == "Chr_X",]

outs_new <- list()
for (i in 1:nrow(data)){
  x <- data[i,]
  outs_new[[i]] <- del(x)
}
outs_new <- outs_new %>% bind_rows()

x <- dif_exp_data[1,]
del <- function(x){
  inv_points <- filter(new_invs2, 
                       start < x[,17] & end > x[,17] | 
                         start < x[,18] & end > x[,18] | 
                         start < x[,17] & end > x[,18] | 
                         start > x[,17] & end < x[,18])
  if (nrow(inv_points) > 0){
    out <- cbind(x, inv_points)
    return(out)
  }
}

del(x)

x <- new_inversions[1,]

ortho_sum <- ortholog_table %>% 
  group_by(REF_GENE_NAME, chr, consensus_gene) %>% 
  summarise(start = min(start), 
            end = max(end)) %>% 
  unique()

filter_invs_func <- function(x){
  a <- filter(ortho_sum, 
         chr == "Chr_X",
         start <= x$start & end >=x$start2 |
           start >= x$start & end <= x$start2 | 
           start <= x$start & end >= x$start | 
           start <=x$start2 & end >= x$start2)
  b <- filter(ortho_sum, 
              chr == "Chr_X",
              start <= x$end & end >=x$end2 |
                start >= x$end & end <= x$end2 | 
                start <= x$end & end >= x$end | 
                start <=x$end2 & end >= x$end2)
  out <- rbind(a,b)
  out$inversion_breakpoint <- x$inversions
  return(out)
}

outs <- list()
for (i in 1:nrow(new_inversions)){
  outs[[i]] <- filter_invs_func(new_inversions[i,])
}
outs <- bind_rows(outs) %>% 
  mutate(REF_GENE_NAME = gsub("_", "-", REF_GENE_NAME))
unique(outs)
outs2 <- outs %>% 
  group_by(REF_GENE_NAME, consensus_gene) %>% 
  summarise(start = min(start), 
            end = max(end)) %>% 
  mutate(search_term = 
           paste0("Chr_X:", start, "-", end)) %>% 
  mutate(REF_GENE_NAME = gsub("_", "-", REF_GENE_NAME))


send_df <- filter(dif_exp_data, genes %in% outs2$REF_GENE_NAME) %>% 
  dplyr::select(genes, consensus_gene, celltype,  Significant, FDR, logFC) %>% 
  rename(Gene = genes, Ortholog = consensus_gene, Differential_expression = Significant, 
         Cell_type = celltype) %>% 
  merge(outs, by.x = c("Gene", "Ortholog"), by.y = c("REF_GENE_NAME", "consensus_gene"), all = T) %>% 
 # dplyr::select(-c(start, end, inversion_breakpoint, chr)) %>% 
  unique()


send_df[is.na(send_df)] <- "Not expressed"
send_df %>% View()

write.table(send_df, "data/breakpoint_alignments/break_point_expression.csv", 
            quote = F, sep = ',', row.names = F)

gois <- intersect(unique(send_df$Gene[send_df$Differential_expression != "Not expressed"]), rownames(seurat_final@assays$RNA$counts))

seurat_final@assays$RNA$counts[gois, ] %>% rowMeans()


