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


load("data/RData/seurat_final.RData")
seurat_final <- JoinLayers(seurat_final)

seurat_final$celltype[seurat_final$celltype %in% c("Early spermatids", "Late spermatids")] <- "Spermatids"
Idents(seurat_final) <- "celltype"

levels(seurat_final) <- c("Muscle", "Early cyst", "Late cyst", 
                          "GSC/Spermatogonia", "Primary spermatocytes", 
                          "Secondary spermatocytes", "Spermatids")
seurat_final@meta.data <- seurat_final@meta.data %>%  
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids")))

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

cell_numbers <- table(seurat_final$celltype, seurat_final$sample)
write.csv(cell_numbers, "data/merged_spermatids_cell_numbers.csv")

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
  cpm <- edgeR::cpm(y, log = F)
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
  out$logcpm <- edgeR::cpm(y, log = T, prior.count = .0001)
  out$cpm <- edgeR::cpm(y, log = F)
  out$toptags <- topTags(res, n = nrow(res))
  out$exacttest <- topTags(STSR, n = nrow(STSR))
  out$func <- func
  out$y <- y
  return(out)
}



de.results_act_edger <- sapply(unique(agg_cell_types$celltype), PBDGE, sce = agg_cell_types, func = 'custom',
                               USE.NAMES = T, simplify = F)
names(de.results_act_edger) <- unique(agg_cell_types$celltype)

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
                                                "Spermatids"))) %>% merge(medianAs) %>% 
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
save(XA, tidy_cpm_edger, dif_exp_data, de.results_act_edger, file = "data/RData/Spermatid_combined_DEG_DC.RData")

dif_exp_data_table <-  dif_exp_data %>% 
  mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  group_by(celltype, Significant, chr) %>% 
  mutate(celltype = gsub("\n", "", celltype)) %>% 
  summarise(count = n()) %>% 
  reshape2::dcast(celltype + chr ~ Significant, value.var = "count") %>% 
  rename(`Cell type` = celltype, Chromosome = chr) %>% 
  .[c(7,8,1,2,5,6,3,4,9,10,11,12,13,14),]

dif_exp_data_table[is.na(dif_exp_data_table)] <- 0


##############################################################

#glm of DGE data #
model1 <- dif_exp_data %>% 
  mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>%
  glm(Significant ~chr, family = binomial, data = .)

model2 <- dif_exp_data %>% 
  mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>%
  glm(Significant ~celltype, family = binomial, data = .)

model3 <- dif_exp_data %>% 
  mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>%
  glm(Significant ~celltype + chr, family = binomial, data = .)


model4 <- dif_exp_data %>% 
  mutate(Significant = ifelse(Significant == "Unbiased", 0, 1)) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst",
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>%
  glm(Significant ~(chr * celltype) , family = binomial, data = .)

summary(model4)

anova(model1, model3, test = "Chisq")
anova(model2, model3, test = "Chisq")
anova(model2, model3, test = "Chisq")
anova(model3, model4, test = "Chisq")

#Model 4 is best fit 
summary(model4)$coefficients %>% 
  write.table("data/DEG_model4_coefficients.tsv", sep = "\t", quote = F, row.names = T)

tmp <- dif_exp_data %>% 
  #filter(Significant != "Unbiased") %>%
  mutate(consensus_gene = ifelse(consensus_gene == genes, " ", consensus_gene)) %>%
  dplyr::select(genes, chr, logFC, FDR, celltype, Significant, consensus_gene) %>% 
  rename("Gene" = genes, "Log2(FC)" = logFC, "FDR" = FDR, 
         "Cell type" = celltype, "Bias" = Significant, 
         "Drosophila ortholog" = consensus_gene) %>% 
  .[,c(1,7,2,5,3,4,6)]



celltypes <- c("M", "EC", "LC", "GSC", "PS", "SS", "ST")
names(celltypes) <- c("Muscle", "Early cyst", "Late cyst", "GSC/Spermatogonia", 
                      "Primary spermatocytes", "Secondary spermatocytes", "Spermatids")

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
sig_dif_exp1[is.na(sig_dif_exp1)] <- "No expression"
keep <- which(rowSums(sig_dif_exp1 == "SR-biased" | sig_dif_exp1 == "ST-biased") > 0)
sig_dif_exp1 <- sig_dif_exp1[keep,]


write.table(sig_dif_exp1, "data/Spermatids_merged_DEG_full_table.tsv", sep = "\t", quote = F, row.names = F)


dif_exp_data_table %>% 
  write.table(., "data/Spermatids_merged_DEG_table.tsv", sep = "\t", quote = F, row.names = F)


## AUTOSOMAL AND X-LINKED GENES ##
Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()

nFeatures_RNA_autosomes <- colSums(seurat_final@assays$RNA$counts[Agenes,] > 0)

seurat_final <- AddMetaData(seurat_final, nFeatures_RNA_autosomes, 'nFeature_RNA_autosomes')

X_exp <- PercentageFeatureSet(seurat_final, features = Xgenes)
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')

Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)/
  colSums(seurat_final@assays$RNA$counts > 1)
Xgene_prop2 <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)/
  nrow(ortholog_table[ortholog_table$chr == "Chr_X",])



seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')
seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop2, col.name = 'Xprop2')
rm_axis = theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.x=element_blank())


######################### MAIN FIGURE PLOTS ######################### 

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7", 'darkgrey')

expected_X_exp <- filter(ortholog_table, chr == "Chr_X") %>%
  dplyr::select(REF_GENE_NAME) %>% 
  unique() %>% nrow()
expected_A_exp <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2")) %>%
  dplyr::select(REF_GENE_NAME) %>% 
  unique() %>% nrow()
expected_expression = expected_X_exp/(expected_A_exp + expected_X_exp)


nfeatures_figure <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA_autosomes, fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + scale_fill_manual(values= cbPalette) + 
  labs(x = "", y = "N\u00b0 expressed  \n autosomal genes") + 
  theme(legend.position="none", 
        axis.line.x = element_blank()) + rm_axis


UMAP <- DimPlot(seurat_final, cols = cbPalette, group.by = 'celltype', label = T) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(legend.position="none")
UMAP_final <- UMAP[[1]]$data %>% 
  ggplot(aes(x = umap_1, y = umap_2, colour = celltype)) + 
  geom_point(stroke = NA, size = .7) + scale_color_manual(values= cbPalette) + 
  guides(colour = guide_legend(override.aes = list(size=3), title = "Cell type")) + 
  theme_classic() + labs(x = "UMAP 1", y = "UMAP 2", legend = 'Cell type') +
  theme(legend.position = 'none', legend.box = "vertical", 
        axis.text.x = element_text(color="black"), 
        axis.ticks = element_line(color = "black")) +
  shadowtext::geom_shadowtext(data = UMAP[[1]]$layers[[2]]$data, aes(label = celltype), colour = 'black', 
                              size = 3, nudge_y = -0.5, nudge_x = 0, 
                              bg.color = 'white', bg.r = 0.1, bg.size = 0.5, 
                              fontface = "bold") + 
  geom_point(data = UMAP[[1]]$layers[[2]]$data, aes(x = umap_1, y = umap_2+0.1), colour = 'black',
             stroke = 1, size = 2, shape = 21, fill = 'white')


# Perform t-test for each group and store results in a data frame
ST_wilcox_results <- XA %>%
  filter(Chr == "X") %>% 
  filter(logcpm >= 1) %>% 
  filter(treatment == "ST") %>% 
  group_by(celltype, treatment) %>%
  summarise(
    p_value = wilcox.test(XA, mu = 0, alternative = 'two.sided')$p.value
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.00001 ~ "***",
      p_value < 0.001 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ " "
    )
  ) %>% 
  mutate(treatment = "ST")


DC_plot <- XA %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = (XA), fill = treatment)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black', linewidth = 0.7) + 
  
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(-5,-1,0.0,5.0), limits = c(-5,8.5)) + 
  guides(fill = 'none') +
  ylab(expression(log[2]~X:A~expression~ratio)) + 
  geom_text(data = ST_wilcox_results, aes(x = celltype, y = 8, label = significance),
            vjust = -0.5, size = 5)  + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black"))

# XCI Plot #

XCI_wilcox_results <- seurat_final@meta.data %>%
  filter(treatment == "ST") %>% 
  group_by(celltype, treatment) %>%
  summarise(
    p_value = wilcox.test(Xprop, mu = expected_expression, alternative = 'two.sided')$p.value
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.00001 ~ "***",
      p_value < 0.001 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ " "
    )
  )



#traj plots 
cbPalette <- c("yellow3", "#0072B2", "#D55E00", "#CC79A7", "darkgrey")


DC_plot_SR <-  XA %>% 
  filter(Chr == "X") %>% 
  filter(logcpm >= 1) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = (XA), fill = treatment)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", fill = "") + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 7, method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
                      symnum.args = list(
                        cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
                        symbols = c("***", "**", "*", " ")), size = 7) + 
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black') + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  theme(legend.position="top", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) +
  ylab(expression(log[2]~(X:A)~expression~ratio)) + 
  scale_y_continuous(breaks=c(-5,0,5)) 


dif_exp_figure <- dif_exp_data %>%
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate(chr = ifelse(chr == "Chr_X", "X chromosome", "Autosomes")) %>%
  mutate(Significant = ifelse(Significant != 'Unbiased', "Biased", "Unbiased")) %>% 
  dplyr::count(celltype, Significant, chr) %>%
  group_by(celltype, chr) %>%
  mutate(Proportion = 100 *(n / sum(n))) %>%
  ungroup() %>% 
  filter(Significant != 'Unbiased') %>% 
  ggplot(aes(x = celltype, y = Proportion, fill = chr)) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_fill_grey(start = 0.3) + 
  theme_classic() + labs(x = "", y = "% of genes differentially expressed\nbetween SR & ST", fill = "") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black"), 
        legend.position = 'top')
########### ALL PLOTS READ ###################





############## MAKING FIGURES ################



A <- UMAP_final +#remove legend
  theme(legend.position = "none")
B <- DC_plot
C <- DC_plot_SR
D <- dif_exp_figure

cts <- seurat_final@meta.data$celltype %>% levels()
comparisons <- list(cts[4:5], cts[5:6], cts[6:7])

E <- seurat_final@meta.data %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  filter(treatment == "ST") %>% 
  ggplot(aes(x = celltype, y = Xprop, fill = treatment)) + geom_violin() +
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  theme(legend.position="none", legend.box = "vertical", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) +
  labs(y = "N\u00b0 expressed X genes / \n N\u00b0 total expressed genes", x = "") + 
  stat_compare_means(comparisons = comparisons, aes(label = ..p.signif..), 
                     label.x = 1.5, label.y = 0.23, method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
                     symnum.args = list(
                       cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
                       symbols = c("***", "**", "*", "ns")), size = 4, 
                     step.increase = 0.04, tip.length = 0.005, ) + 
  ylim(0,0.29) 
Fplot <- seurat_final@meta.data %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  ggplot(aes(x = celltype, y = Xprop, fill = treatment)) + geom_boxplot(outlier.alpha = 0.1) +
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  theme(legend.position="right", legend.box = "vertical", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) +
  labs(y = "N\u00b0 expressed X genes / \n N\u00b0 total expressed genes", x = "", fill = '') + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 0.25, method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
                      symnum.args = list(
                        cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
                        symbols = c("***", "**", "*", " ")), size = 4)

  
  
part1 <- ggarrange(A,D,B,E, labels = c("A", "B", "C", "D"))
part2 <- ggarrange(C, Fplot, common.legend = TRUE,
                   labels = c("E", "F"), legend = 'top')

Extra_plot <- ggarrange(part1, part2, nrow = 2, heights = c(2,1))
ggsave("plots/manuscript_plots/SX_merged_spermatids.pdf", 
       plot = Extra_plot, 
       width = 10, height = 16, units = "in", dpi = 300)
system("open plots/manuscript_plots/SX_merged_spermatids.pdf")


Extra_plot <- ggarrange(A,D,B,C,E,Fplot,
          ncol = 2, nrow = 3,
          common.legend = TRUE, legend = "right",
          labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 20, face = "bold"))











