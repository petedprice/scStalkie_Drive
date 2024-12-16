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
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()

seurat_final <- JoinLayers(seurat_final)
sce <- seurat_final %>% as.SingleCellExperiment(assay = "RNA")

################################################################################


######################## PSEUDOBULK INDIVIDUAL CELL TYPES ######################
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]


## NEW PBDGE FUNCTION
PBDGE <- function(label, sce, func = 'custom', seurat_obj = seurat_final, consistent_genes = T){
  if (consistent_genes == TRUE){
    ST_cells <- filter(seurat_final@meta.data, celltype == label & treatment == "ST") %>% 
      dplyr::select(cells) %>% c() %>% unlist()
    SR_cells <- filter(seurat_final@meta.data, celltype == label & treatment == "SR") %>% 
      dplyr::select(cells) %>% c() %>% unlist()
    ST_keep <- (rowSums(seurat_final@assays$RNA$counts[,ST_cells] > 1) > (0.05 * length(ST_cells))) #| 
                  #rowSums(seurat_final@assays$RNA$counts[,ST_cells] > 1) > 10)
    
    SR_keep <- (rowSums(seurat_final@assays$RNA$counts[,SR_cells] > 1) > (0.05 * length(SR_cells)))# | 
#                  rowSums(seurat_final@assays$RNA$counts[,SR_cells] > 1) > 10)
    
    keep <- ST_keep | SR_keep
  }else {
    keep <- rep(T, nrow(sce))
  }
  dge_data <- sce[keep,sce$celltype == label]
  y <- DGEList(counts(dge_data), samples=colData(dge_data), group = dge_data$treatment)
  y <- calcNormFactors(y)
  design <- model.matrix(~treatment, y$samples)
  y <- normLibSizes(y)
  y <- estimateDisp(y, design)
  cpm <- cpm(y, log = F)
  if (func == "custom"){
    keep1 <- cpm[,y$samples$treatment== "ST"] >= 1
    keep1 <- rowSums(keep1) > (ncol(keep1)/2)
    keep2 <- cpm[,y$samples$treatment== "SR"] >= 1
    keep2 <- rowSums(keep2) > (ncol(keep2)/2)
    keep <- keep1 | keep2
  } else if (func == "edger"){  
    keep <- filterByExpr(y, group = y$sample$treatment, 
                         large.n = 3, min.prop = 0.6, min.count = 5)
    #keep <- filterByExpr(y, group = y$sample$treatment)
  }
  y <- y[keep,keep.lib.sizes=FALSE]
  
  STSR <- exactTest(y, pair = c("ST", "SR")) 
  fit <- glmQLFit(y, design, robust = T)
  res <- glmQLFTest(fit, coef=ncol(design))
  out <- list()
  out$res <- res
  out$logcpm <- cpm(y, log = T, prior.count = .0001)
  out$cpm <- cpm(y, log = F)
  out$toptags <- topTags(res, n = nrow(res))
  out$exacttest <- topTags(STSR, n = nrow(STSR))
  out$func <- func
  out$y <- y
  return(out)
}



de.results_act_edger <- sapply(unique(agg_cell_types$celltype), PBDGE, sce = agg_cell_types, func = 'custom',
                               USE.NAMES = T, simplify = F)


cols <- c("treatment", "celltype", "sample")
get_cpm_values_func <- function(x){
  md <- x$res
  ST <- rownames(filter(md$samples, treatment == "ST"))
  SR <- rownames(filter(md$samples, treatment == "SR"))
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
tidy_cpm_edger$func <- "custom"

################ Medians approach for dosage compensation ------
medianAs <- tidy_cpm_edger %>% 
  group_by(celltype, treatment, genes, Chr) %>%
  summarise(logcpm = log2(mean(2^logcpm))) %>% 
  group_by(celltype, treatment) %>%
  summarise(medianA = median(logcpm[Chr == "A" & logcpm >= 1]), 
            n = length(logcpm[Chr == "A" & logcpm >= 1]))


XA <- tidy_cpm_edger %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% merge(medianAs) %>% 
  group_by(celltype, genes, Chr, treatment, medianA) %>% 
  summarise(logcpm = log2(mean(2^logcpm))) %>% 
  #group_by(sample, celltype, treatment, func, Chr) %>%
  summarise(XA = logcpm - medianA, 
            logcpm = logcpm )



## No of DGE ##


########################  ENRICHMENT FIGURE ##########################
get_DGE_data <- function(ct,x){out = x[[ct]][['toptags']]$table
                             out$celltype <-  ct
                             out$genes <- rownames(out)
                             rownames(out) <- NULL
                             return(out)}

dif_exp_data <- lapply(names(de.results_act_edger), get_DGE_data, x = de.results_act_edger) %>% 
  bind_rows() %>% 
  mutate(Significant = ifelse(.$logFC > 1 & .$FDR < 0.05, "ST-biased",
                                                              ifelse(.$logFC < -1 & .$FDR < 0.05, "SR-biased", "Unbiased"))) %>% 
  merge(ortholog_table, by.x = 'genes', by.y = 'REF_GENE_NAME')
save(XA, tidy_cpm_edger, dif_exp_data, de.results_act_edger, file = "data/RData/DEG_DC.RData")

dif_exp_data_table <-  dif_exp_data %>% 
  mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  group_by(celltype, Significant, chr) %>% 
  mutate(celltype = gsub("\n", "", celltype)) %>% 
  summarise(count = n()) %>% 
  reshape2::dcast(celltype + chr ~ Significant, value.var = "count") %>% 
  rename(`Cell type` = celltype, Chromosome = chr) %>% 
  .[c(3,4,7,8,5,6,1,2,9:14),]

dif_exp_data_table[is.na(dif_exp_data_table)] <- 0


##############################################################

#glm of DGE data #
model1 <- dif_exp_data %>% 
  mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  #mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>%
  glm(Significant ~chr, family = binomial, data = .)

model2 <- dif_exp_data %>% 
  mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  #mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>%
  glm(Significant ~1 +celltype, family = binomial, data = .)

model3 <- dif_exp_data%>% 
  mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  #mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>%
  glm(Significant ~celltype + chr, family = binomial, data = .)


model4 <- dif_exp_data %>% 
 mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  #mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>%
  glm(Significant ~(chr * celltype), family = binomial, data = .)

anova(model1, model3, test = "Chisq")
anova(model2, model3, test = "Chisq")
anova(model2, model3, test = "Chisq")
anova(model3, model4, test = "Chisq")

#Model 3 is best fit 
summary(model3)$coefficients %>% 
  write.table("data/DEG_model3_coefficients.tsv", sep = "\t", quote = F, row.names = T)

tmp <- dif_exp_data %>% 
  #filter(Significant != "Unbiased") %>%
  mutate(consensus_gene = ifelse(consensus_gene == genes, " ", consensus_gene)) %>%
  dplyr::select(genes, chr, logFC, FDR, celltype, Significant, consensus_gene) %>% 
  rename("Gene" = genes, "Log2(FC)" = logFC, "FDR" = FDR, 
         "Cell type" = celltype, "Bias" = Significant, 
         "Drosophila ortholog" = consensus_gene) %>% 
  .[,c(1,7,2,5,3,4,6)]



celltypes <- c("M", "EC", "LC", "GSC", "PS", "SS", "EST", "LST")
names(celltypes) <- c("Muscle", "Early cyst", "Late cyst", "GSC/Spermatogonia", 
                                   "Primary spermatocytes", "Secondary spermatocytes", "Early spermatids", "Late spermatids")

make_df_func <- function(tmp, celltype){
  tmp2 <- tmp %>% 
    filter(`Cell type` == celltype) %>% 
    .[,c(1,2,3,5,6,7)]
  colnames(tmp2)[c(4,5,6)] <- paste(celltypes[celltype], colnames(tmp2)[c(4,5,6)], sep = " ")
  tmp2[,c(4,5)] <- round(tmp2[,c(4,5)], 3)
  return(tmp2)
}

sig_dif_exp1 <- lapply(names(celltypes), make_df_func, tmp = tmp) %>% 
  purrr::reduce(full_join, by = c("Gene", "chr", "Drosophila ortholog"))
sig_dif_exp1[is.na(sig_dif_exp1)] <- "Low Exp"
keep <- which(rowSums(sig_dif_exp1 == "SR-biased" | sig_dif_exp1 == "ST-biased") > 0)
sig_dif_exp1 <- sig_dif_exp1[keep,]


write.table(sig_dif_exp1, "data/DEG_full_table.tsv", sep = "\t", quote = F, row.names = F)


dif_exp_data_table %>% 
  .[c(3,4,7,8,5,6,1,2,9:14),] %>% 
  write.table(., "data/DEG_table.tsv", sep = "\t", quote = F, row.names = F)







################# INVERSION CHECK --------------
DGE_inversion <- dif_exp_data

inversions <- data.frame(inversions = c("inv1", "inv2", "inv3", "inv4", "inv5"),
                         starts = c(76.7, 80.2, 59,32,1.2), 
                         ends = c(80.4, 87.07, 76.7, 46.3, 9.03)) %>% 
  mutate(starts = starts * 1e6, 
         ends = ends * 1e6)

inv_func <- function(x, inversions){
  if (x[12] == "Chr_X"){
    st <- inversions[,1][which(inversions[,2] < as.numeric(x[16]) & inversions[,3] > 
                                 as.numeric(x[17]))]
    if(length(st) == 0){
      return("outside inversion")
    } else {
      return(st[1])
      #  return("inversion")
    }
  }else {
    return("autosomal")
  }
}


DGE_inversion$inversion <- apply(DGE_inversion, 1, inv_func, inversions = inversions) %>% unlist()
chisq_inv_table <- DGE_inversion %>% filter(chr == "Chr_X") %>% 
  group_by(inversion, Significant, celltype) %>% 
  summarise(count = n()) %>% 
  reshape2::dcast(celltype + inversion ~ Significant, value.var = "count") 
chisq_inv_table[is.na(chisq_inv_table)] <- 0


subset_chisq <- function(cell, data){
  ss <- filter(data, celltype == cell) %>% 
    dplyr::select(-c(celltype, inversion)) %>% 
    chisq.test()
  
  return(ss)
}

inv_chisqs <- lapply(unique(chisq_inv_table$celltype), subset_chisq, data = chisq_inv_table)

model_inv <- DGE_inversion %>% 
  filter(chr == "Chr_X") %>% 
  mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  #mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>%
  glm(Significant ~celltype + inversion, family = binomial, data = .)

summary(model_inv)$coefficients %>% 
  write.table("data/DEG_model_inv_coefficients.tsv", sep = "\t", quote = F, row.names = T)
#There are no interesting results for enrichment of DGE inside/outside the inversion 
