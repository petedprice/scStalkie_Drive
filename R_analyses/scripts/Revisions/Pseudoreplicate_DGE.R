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
library(lme4)
library(GenomicFeatures)
library(ggrepel)
library(lme4)
rm(list = ls())

load("data/RData/seurat_final.RData")
DefaultAssay(seurat_final) <- "RNA"
seurat_final$treatment <- "SR"
seurat_final$treatment[grep("st", seurat_final$sample)] <- "ST"
ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene, FBconcensus) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()

seurat_final <- JoinLayers(seurat_final)

################################################################################


######################## PSEUDOBULK INDIVIDUAL CELL TYPES ######################
#Add subsample
seurat_final$subsample <- NA
for (s in unique(seurat_final$sample)){
  for (ct in unique(seurat_final$celltype)){
    no_cells <- nrow(filter(seurat_final@meta.data, celltype == ct, sample == s))
    subsamples <- paste0(s, "_", sample(c("a", "b", "c", "d"), no_cells, replace = T))
    seurat_final$subsample[seurat_final$celltype == ct & seurat_final$sample == s] <- 
      subsamples
    
  }
}
sce <- seurat_final %>% as.SingleCellExperiment(assay = "RNA")



agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "subsample")])
table(agg_cell_types$celltype, agg_cell_types$subsample)
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]
table(agg_cell_types$celltype, agg_cell_types$sample)


celltype = "GSC/Spermatogonia"
s1 = 'sr1'
s2 = "sr2"
func = 'custom'
sce = agg_cell_types
consistent_genes = T
func = 'custom'
## NEW PBDGE FUNCTION
PBDGE <- function(celltype, s1, s2, sce, func = 'custom', seurat_obj = seurat_final, consistent_genes = T){
  if (consistent_genes == TRUE){
    group1_cells <- filter(seurat_final@meta.data, sample == s1 & celltype == ct) %>% 
      dplyr::select(cells) %>% c() %>% unlist()
    group2_cells <- filter(seurat_final@meta.data, sample == s2 & celltype == ct) %>% 
      dplyr::select(cells) %>% c() %>% unlist()
    g1_keep <- (rowSums(seurat_final@assays$RNA$counts[,group1_cells] > 1) > (0.05 * length(group1_cells))) #| 
    #rowSums(seurat_final@assays$RNA$counts[,ST_cells] > 1) > 10)
    
    g2_keep <- (rowSums(seurat_final@assays$RNA$counts[,group2_cells] > 1) > (0.05 * length(group2_cells)))# | 
    #                  rowSums(seurat_final@assays$RNA$counts[,SR_cells] > 1) > 10)
    
    keep <- g1_keep | g2_keep
  }else {
    keep <- rep(T, nrow(sce))
  }
  dge_data <- sce[keep,sce$sample %in% c(s1, s2) & sce$celltype == celltype]
  y <- DGEList(counts(dge_data), samples=colData(dge_data), group = dge_data$sample)
  y <- calcNormFactors(y)
  design <- model.matrix(~sample, y$samples)
  y <- normLibSizes(y)
  y <- estimateDisp(y, design)
  cpm <- edgeR::cpm(y, log = F)
  if (func == "custom"){
    keep1 <- cpm[,y$samples$sample== s1] >= 1
    keep1 <- rowSums(keep1) > (ncol(keep1)/2)
    keep2 <- cpm[,y$samples$sample== s2] >= 1
    keep2 <- rowSums(keep2) > (ncol(keep2)/2)
    keep <- keep1 | keep2
  } else if (func == "edger"){  
    keep <- filterByExpr(y, group = y$sample$treatment, 
                         large.n = 3, min.prop = 0.6, min.count = 5)
    #keep <- filterByExpr(y, group = y$sample$treatment)
  }
  y <- y[keep,keep.lib.sizes=FALSE]
  
  s1s2 <- exactTest(y, pair = c(s1, s2)) 
  fit <- glmQLFit(y, design, robust = T)
  res <- glmQLFTest(fit, coef=ncol(design))
  out <- list()
  out$res <- res
  out$logcpm <- edgeR::cpm(y, log = T, prior.count = .0001)
  out$cpm <- edgeR::cpm(y, log = F)
  out$toptags <- topTags(res, n = nrow(res))
  out$exacttest <- topTags(s1s2, n = nrow(s1s2))
  out$func <- func
  out$y <- y
  return(out)
}

de.results_act_edger <- list()

run <- F
if (run == T){
  for (ct in unique(agg_cell_types$celltype)){
    samples <- unique(agg_cell_types$sample)
    
    print(paste0("running all comparisons for celltype ", ct))
    for (s1 in samples){
      print(paste0("running all comparisons for ", s1))
      for (s2 in samples){
        if (s1 == s2){
          next
        } else {
          print(s2)
          if (sum(agg_cell_types$sample == s1 & agg_cell_types$celltype == ct) < 2| 
              sum(agg_cell_types$sample == s2 & agg_cell_types$celltype == ct) < 2) {
            print('not enough samples for this comparison')
            next
          } else {
            tmp <- PBDGE(ct, s1, s2, agg_cell_types)
            de.results_act_edger[[paste(ct, s1, s2, sep = "_")]] <- tmp
          }
        }
        samples <- samples[samples != s1]
      }
    }
  }
  
  
  
  get_DGE_data <- function(ct,s1,s2,x){
    out = x[[paste(ct, s1, s2, sep = "_")]][['toptags']]$table
    out$celltype <-  ct
    out$s1 <- s1
    out$s2 <- s2
    out$genes <- rownames(out)
    rownames(out) <- NULL
    ST <- grepl("st", c(s1,s2))
    SR <- grepl("sr", c(s1,s2))
    if (sum(ST) == 1){
      comp <- "STvsSR"
    } else if (sum(ST) == 2){
      comp <- "STvsST"
    } else {
      comp = "SRvsSR"
    }
    out$comparison <- comp   
    return(out)
  }
  
  samples <- unique(agg_cell_types$sample)
  dif_exp_data_list <- list()
  for (ct in unique(agg_cell_types$celltype)){
    samples <- unique(agg_cell_types$sample)
    print(paste0("running all comparisons for celltype ", ct))
    for (s1 in samples){
      print(paste0("running all comparisons for ", s1))
      for (s2 in samples){
        if (s1 == s2){
          next
        } else {
          print(s2)
          tmp <- get_DGE_data(ct, s1, s2, de.results_act_edger)
          dif_exp_data_list[[paste(ct, s1, s2, sep = "_")]] <- tmp
        }
        samples <- samples[samples != s1]
      }
    }
  }
  
  
  dif_exp_data <- dif_exp_data_list %>% 
    bind_rows() %>% 
    mutate(Significant = ifelse(.$logFC > 1 & .$FDR < 0.05, "s1-biased",
                                ifelse(.$logFC < -1 & .$FDR < 0.05, "s2-biased", "Unbiased"))) %>% 
    merge(ortholog_table, by.x = 'genes', by.y = 'REF_GENE_NAME')
  save(dif_exp_data, file = 'data/pseudo_replicate.RData')
  pr_dif_exp_data <- dif_exp_data
} else {
  load("data/pseudo_replicate.RData")
  pr_dif_exp_data <- dif_exp_data
}  

load("data/RData/DEG_DC.RData")

A <- pr_dif_exp_data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate(comp2 = paste0(s1,s2)) %>% 
  filter(comparison == "STvsST") %>% 
  mutate(chr = ifelse(chr == "Chr_X", "X chromosome", "Autosomes")) %>%
  mutate(Significant = ifelse(Significant != 'Unbiased', "Biased", "Unbiased")) %>% 
  dplyr::count(celltype, Significant, chr, comp2) %>%
  group_by(celltype, chr, comp2) %>%
  mutate(Proportion = 100 *(n / sum(n))) %>%
  ungroup() %>% 
  filter(Significant != 'Unbiased') %>% 
  ggplot(aes(x = comp2, y = n)) +geom_boxplot()+ ylim(0,4000)
  # ggplot(aes(x = celltype, y = Proportion, fill = chr)) + 
  # geom_boxplot() + 
  # scale_fill_grey(start = 0.3) + 
  # theme_classic() + labs(x = "", y = "% of genes differentially expressed\nbetween SR & ST", fill = "") + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
  #       axis.ticks = element_line(color = "black"), 
  #       legend.position = 'top') + 
  # labs(title = "STvsST")


B <- pr_dif_exp_data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate(comp2 = paste0(s1,s2)) %>% 
  filter(comparison == "SRvsSR") %>% 
  mutate(chr = ifelse(chr == "Chr_X", "X chromosome", "Autosomes")) %>%
  mutate(Significant = ifelse(Significant != 'Unbiased', "Biased", "Unbiased")) %>% 
  dplyr::count(celltype, Significant, chr, comp2) %>%
  group_by(celltype, chr, comp2) %>%
  mutate(Proportion = 100 *(n / sum(n))) %>%
  ungroup() %>% 
  filter(Significant != 'Unbiased') %>% 
  ggplot(aes(x = comp2, y = n)) + geom_boxplot()+ ylim(0,4000)
  # ggplot(aes(x = celltype, y = Proportion, fill = chr)) + 
  # geom_boxplot() + 
  # scale_fill_grey(start = 0.3) + 
  # theme_classic() + labs(x = "", y = "% of genes differentially expressed\nbetween SR & ST", fill = "") + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
  #       axis.ticks = element_line(color = "black"), 
  #       legend.position = 'top') + 
  # labs(title = "SRvsSR")

C <- pr_dif_exp_data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate(comp2 = paste0(s1,s2)) %>% 
  filter(comparison == "STvsSR") %>% 
  mutate(chr = ifelse(chr == "Chr_X", "X chromosome", "Autosomes")) %>%
  mutate(Significant = ifelse(Significant != 'Unbiased', "Biased", "Unbiased")) %>% 
  dplyr::count(celltype, Significant, chr, comp2) %>%
  group_by(celltype, chr, comp2) %>%
  mutate(Proportion = 100 *(n / sum(n))) %>%
  ungroup() %>% 
  filter(Significant != 'Unbiased') %>% 
  ggplot(aes(x = comp2, y = n)) + geom_boxplot() + ylim(0,4000)
  # ggplot(aes(x = celltype, y = Proportion, fill = chr)) + 
  # geom_boxplot() + 
  # scale_fill_grey(start = 0.3) + 
  # theme_classic() + labs(x = "", y = "% of genes differentially expressed\nbetween SR & ST", fill = "") + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
  #       axis.ticks = element_line(color = "black"), 
  #       legend.position = 'top') + 
  # labs(title = "SRvsST")

ggarrange(A,B, C, ncol = 3)  










plot <- pr_dif_exp_data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate(comp2 = paste0(s1,s2)) %>% 
  filter(comparison != "SRvsSR")%>% 
  mutate(chr = ifelse(chr == "Chr_X", "X chromosome", "Autosomes")) %>%
  mutate(Significant = ifelse(Significant != 'Unbiased', "Biased", "Unbiased")) %>% 
  dplyr::count(comparison, celltype, Significant, chr, comp2) %>%
  group_by(celltype, chr, comp2, comparison) %>%
  mutate(Proportion = 100 *(n / sum(n))) %>%
  ungroup() %>% 
  filter(Significant != 'Unbiased') %>% 
  ggplot(aes(x = comparison, y = Proportion, fill = chr)) + 
  geom_boxplot() + 
  scale_fill_grey(start = 0.4, end = 0.7) + 
  labs(fill = "", y = "% of genes differentially expressed\nbetween SR & ST") + 
  theme_classic() + 
  facet_wrap(~celltype, ncol = 2)
#  scale_fill_brewer(palette = "Grayscale")
plot





