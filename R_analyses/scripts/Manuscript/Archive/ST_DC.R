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

load("data/RData/ST_seurat_final.RData")
DefaultAssay(ST_seurat_final) <- "RNA"
ST_seurat_final$treatment <- "SR"
ST_seurat_final$treatment[grep("st", ST_seurat_final$sample)] <- "ST"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(ST_seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(ST_seurat_final))%>% unique()

sce <- ST_seurat_final %>% as.SingleCellExperiment(assay = "RNA")

################################################################################


######################## PSEUDOBULK INDIVIDUAL CELL TYPES ######################
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]


### WORKIGN IN PROGERESS ###


## NEW PBDGE FUNCTION
label = 'Spermatids'
PBDGE <- function(label, sce, func = 'custom', seurat_obj = ST_seurat_final, consistent_genes = T){
  if (consistent_genes == TRUE){
    ST_cells <- filter(ST_seurat_final@meta.data, celltype == label & treatment == "ST") %>% 
      dplyr::select(cells) %>% c() %>% unlist()
    ST_keep <- (rowSums(ST_seurat_final@assays$RNA$counts[,ST_cells] > 1) > (0.05 * length(ST_cells)))

    keep <- ST_keep
  }else {
    keep <- rep(T, nrow(sce))
  }
  dge_data <- sce[keep,sce$celltype == label]
  y <- DGEList(counts(dge_data), samples=colData(dge_data))
  y <- normLibSizes(y)
  y <- estimateDisp(y)
  cpm <- cpm(y, log = F)
  if (func == "custom"){
    keep1 <- cpm[,y$samples$treatment== "ST"] >= 2
    keep1 <- rowSums(keep1) > (ncol(keep1)/2)
    keep <- keep1
  } else if (func == "edger"){  
    keep <- filterByExpr(y, group = y$sample$treatment, 
                         large.n = 3, min.prop = 0.6, min.count = 5)
  }
  y <- y[keep,keep.lib.sizes=FALSE]
  

  out <- list()
  out$y <- y
  out$logcpm <- cpm(y, log = T, prior.count = .0001)
  out$cpm <- cpm(y, log = F)
  out$func <- func
  return(out)
}


de.results_act_edger <- sapply(unique(agg_cell_types$celltype), PBDGE, sce = agg_cell_types, func = 'custom',
                               USE.NAMES = T, simplify = F)


cols <- c("treatment", "celltype", "sample")
get_cpm_values_func <- function(x){
  md <- x$y
  ST <- rownames(filter(md$samples, treatment == "ST"))
  cpm <- x$logcpm %>% 
    as.data.frame() %>%
    rownames_to_column("genes") %>% 
    pivot_longer(cols = colnames(.)[2:ncol(.)], names_to = "sample", values_to = "logcpm") %>% 
    mutate(Chr = ifelse(genes %in% Xgenes, "X", 
                        ifelse(genes %in% Agenes, "A", "Mito"))) %>% 
    mutate(treatment = ifelse(sample %in% ST, "ST", "SR")) %>% 
    mutate(celltype = md$samples$celltype[1]) %>% 
    rename(edgeRsample = sample) %>% 
    mutate(sample = md$samples[edgeRsample,]$sample)
  return(cpm)
}  

tidy_cpm_edger <- lapply(de.results_act_edger, get_cpm_values_func) %>% bind_rows()
tidy_cpm_edger$func <- "edgeR"

################ Medians approach for dosage compensation ------
medianAs <- tidy_cpm_edger %>% 
  group_by(celltype, treatment, genes, Chr) %>%
  summarise(logcpm = log2(mean(2^logcpm))) %>% 
  group_by(celltype, treatment) %>%
  summarise(medianA = median(logcpm[Chr == "A" & logcpm >= 1]), 
            n = length(logcpm[Chr == "A" & logcpm >= 1]))


XA <- tidy_cpm_edger %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% merge(medianAs) %>% 
  group_by(celltype, genes, Chr, treatment, medianA) %>% 
  summarise(logcpm = log2(mean(2^logcpm))) %>% 
  #group_by(sample, celltype, treatment, func, Chr) %>%
  summarise(XA = logcpm - medianA, 
            logcpm = logcpm )

ST_XA <- XA
ST_tidy_cpm_edger <- tidy_cpm_edger

save(ST_XA, ST_tidy_cpm_edger, file = "data/RData/ST_DEG_DC.RData")




