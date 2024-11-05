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
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

load("data/RData/seurat_final.RData")
load("data/RData/DEG_DC.RData")

DefaultAssay(seurat_final) <- "RNA"
seurat_final <- subset(seurat_final, celltype %in% c("GSC/Spermatogonia", "Primary spermatocytes", 
                                                     "Secondary spermatocytes", "Spermatids"))
                                                    

seurat_final@meta.data$celltype_broad <- seurat_final$celltype
seurat_final@meta.data$celltype_broad[seurat_final$celltype %in% c("Secondary spermatocytes", "Spermatids")] <- 
  "Hap"
seurat_final@meta.data$celltype_broad[seurat_final$celltype %in% c("Primary spermatocytes", "GSC/Spermatogonia")] <- 
  "Dip"


sce <- seurat_final %>% as.SingleCellExperiment(assay = "RNA") 
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype_broad", "sample", "treatment")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]
###############################################





############# CREATING DEG DATA ---------------

y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))

y <- calcNormFactors(y)
y <- normLibSizes(y)
y <- estimateDisp(y)
cpm <- cpm(y, log = F)
keep <- rep(F, nrow(cpm))
for (ct in unique(agg_cell_types$celltype_broad)){
  keep1 <- cpm[,y$samples$celltype_broad== ct] >= 2
  keep2 <- rowSums(keep1) >= (ncol(keep1)/2)
  keep <- keep | keep2
}
sum(keep)

keep <- filterByExpr(y, group = y$sample$celltype_broad, 
                      large.n = 3, min.prop = 0.6, min.count = 5)

y <- y[keep,keep.lib.sizes=FALSE]
cpm <- cpm(y, log = F)

#Filter low cell count expressions 
tau_keep_genes <- list()
for (ct in unique(seurat_final$celltype_broad)){
  counts <- seurat_final@assays$RNA$counts
  ST_cells <- filter(seurat_final@meta.data, celltype_broad == ct & 
                       treatment == "ST") %>% 
      dplyr::select(cells) %>% c() %>% unlist()
  ST_keep <- rowSums(counts[,ST_cells] > 1) > (0.05 * length(ST_cells))
  ST_keep_genes <- names(ST_keep)[ST_keep]
  
  SR_cells <- filter(seurat_final@meta.data, celltype_broad == ct & 
    treatment == "SR") %>% 
    dplyr::select(cells) %>% c() %>% unlist()
  SR_keep <- rowSums(counts[,SR_cells] > 1) > (0.05 * length(SR_cells))
  SR_keep_genes <- names(SR_keep)[SR_keep]
  tau_keep_genes[[ct]] <-  intersect(SR_keep_genes, ST_keep_genes)
}


#Calculating Tau
tau <- function(gene, cpm, y){
  expression <- cpm[gene,]
  #expression[expression < 1] <- 0
  #expression <- log2(expression + 0.000001)
  sample <- y$samples$sample
  celltype_broad <- y$samples$celltype_broad
  taudata <- data.frame(exp = expression, sample = sample, 
                        celltype_broad = celltype_broad) %>% 
    group_by(celltype_broad) %>% 
    summarise(mean = mean(exp), 
              log2mean = log2(mean))
  
  taudata$Xhat <- taudata$log2mean/max(taudata$log2mean)
  taudata$tau <- sum(1-taudata$Xhat)/(nrow(taudata)-1)
  taudata$gene <- gene
  taudata$ct_specificity = taudata$celltype_broad[which.max(taudata$Xhat)]
  return(taudata)
  
}

taus <- do.call(rbind, lapply(rownames(cpm), tau, cpm = cpm, y = y))
taus <- as.data.frame(taus)
taus <- taus %>% 
  merge(ortholog_table, by.x = 'gene', by.y = 'REF_GENE_NAME')

tau_sum_func <- function(chrom, ct, taus, out = 'datas', tau_keep_genes = tau_keep_genes){
  tau_genes <- filter(taus, tau >= 0.5)%>% 
    filter(chr == chrom) %>% 
    filter(ct_specificity == ct & 
             celltype_broad == ct) %>% 
    filter(gene %in% tau_keep_genes[[ct]]) %>% 
    dplyr::select(gene, mean) %>% unique()
  tgs <- tau_genes$gene %>% unique()
  tau_exp <- PercentageFeatureSet(seurat_final, features = tgs)
  tau_metadata <- AddMetaData(object = seurat_final, metadata = tau_exp, col.name = 'tau_exp')@meta.data %>% 
    dplyr::select(celltype, celltype_broad, tau_exp, treatment, sample)
  tau_metadata$chr <- chrom
  tau_metadata$tau_type <- ct
  if (out == 'datas'){
    return(tau_metadata)
  } else {
    if (nrow(tau_genes) > 0){
      out_genes <- cbind(tau_genes, chr = chrom, tau_type = ct)
      return(out_genes)
    } else {
      return(NULL)
    }
  }
}

#using the celltype_broads and unique chromosomes do a double lapply and bind_rows of the above function 
chroms <- unique(taus$chr)
cts <- unique(taus$ct_specificity)
tau_data <- lapply(chroms, function(x) lapply(cts, tau_sum_func, chrom = x, taus = taus, tau_keep_genes = tau_keep_genes)) %>% 
  bind_rows()
tau_genes <- lapply(chroms, function(x) lapply(cts, tau_sum_func, chrom = x, taus = taus, out = 'genes', tau_keep_genes = tau_keep_genes)) %>% 
  bind_rows()
table(tau_genes$chr, tau_genes$tau_type)

#save(tau_data, tau_genes, file = "data/RData/tau_data.RData")

rm_axis = theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())
 
tau_plot_data <- tau_data #%>% 
  #mutate(tau_type = paste0(tau_type, " specific genes")) #%>% 
  #mutate(tau_type = factor(tau_type, levels = c("Muscle specific genes", "Cyst specific genes", 
  #                                              "Haploid-germ specific genes", "Diploid-germ specific genes"))) %>% 
  #mutate(chr = gsub("Chr_", "Chromosome  ", chr))



tau_plot_data %>% 
  ggplot(aes(x = celltype, y = tau_exp, fill = treatment)) + geom_boxplot(outlier.alpha = 0.1) +
  labs(y = "Percentage of transcripts from \nPre-meiosis specific genes", x = "") + 
  scale_fill_brewer(palette = 'Set2') +
  facet_wrap(~chr + tau_type, nrow = 3)
  


p1 <- tau_plot_data %>% filter(tau_type == "Dip") %>%
  ggplot(aes(x = celltype, y = tau_exp, fill = treatment)) + geom_boxplot(outlier.alpha = 0.1) +
  labs(y = "Percentage of transcripts from \nPre-meiosis specific genes", x = "") + 
  scale_fill_brewer(palette = 'Set2') +
  facet_wrap(~chr, nrow = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme_classic() + 
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  rm_axis
    #ylim(0,0.8) 
p2 <- tau_plot_data %>% filter(tau_type == "Hap") %>%
  ggplot(aes(x = celltype, y = tau_exp, fill = treatment)) + geom_boxplot(outlier.alpha = 0.1) +
  labs(y = "Percentage of transcripts from \nPost-meiosis specific genes", x = "") + 
  scale_fill_brewer(palette = 'Set2') +
  facet_wrap(~chr, nrow = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme_classic() + 
  theme(legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) #+ 
 # ylim(0,0.8)
  
ggarrange(plotlist = list(p1,p2), nrow = 2, common.legend = TRUE, heights = c(1,1.3))



### FINDALLMARKERS APPROACH ###
Idents(seurat_final) <- "celltype"
markers <- seurat_final %>% 
  #subset(., treatment == "ST")  %>%
  FindMarkers(., ident.1 = 'Hap', ident.2 = 'Dip', 
              only.pos = T)

gois <- markers %>% 
  #filter(cluster == "Dip") %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 2 & pct.1 > 0.1 & 
           pct.2 < 0.1) %>% 
  mutate(gene = row.names(.)) %>% 
  merge(ortholog_table, by.x = 'gene', by.y = 'REF_GENE_NAME') %>% 
  filter(chr == "Chr_X") %>% 
  dplyr::select(gene) %>% unlist()
length(gois)


tau_exp <- PercentageFeatureSet(seurat_final, features = gois)
tau_metadata <- AddMetaData(object = seurat_final, metadata = tau_exp, col.name = 'tau_exp')@meta.data %>% 
  dplyr::select(celltype, celltype_broad, tau_exp, treatment, sample)

tau_metadata %>% 
  ggplot(aes(x = celltype, y = tau_exp, fill = treatment)) + 
  geom_boxplot() + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 1, method = 'wilcox.test', 
                      method.args = list(alternative = 'two.sided'))



#
filter(dif_exp_data, genes %in% gois) %>% 
  filter(celltype == "Spermatids") %>% 
  View()







### Tau DGE ##
tau_DGE <- taus %>% 
  dplyr::select(gene, celltype_broad, mean, Xhat, tau, ct_specificity) %>% 
  merge(dif_exp_data, by.x = 'gene', by.y = 'genes')

tau_DGE %>% 
  filter(celltype_broad == "Post-meiotic") %>%
  filter(ct_specificity == "Post-meiotic") %>% 
  filter(abs(logFC > 1) & FDR < 0.05) %>% 
  filter(celltype %in% c("Secondary spermatocytes", "Spermatids")) %>% 
  View()





