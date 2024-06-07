## based off http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
# load libraries and functions
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


### Data setup ---
load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned_ss.RData")
Muscle <- c(17)
Spermatocytes <- c(10,11,13)
Spermatids <- c(3,4)
`GSC/Spermatogonia` <- c(0,2,7)
Pre_meiotic_cyst <- c(5,12,14,15)
Post_meiotic_cyst <- c(1,6,8,9,16)

seurat_integrated_ss@meta.data$celltype <- 'NA'
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Muscle] <- "Muscle"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatocytes] <- "Spermatocytes"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatids] <- "Spermatids"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% `GSC/Spermatogonia`] <- "GSC & Spermatogonia"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Pre_meiotic_cyst] <- "Pre-meiotic \ncyst"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Post_meiotic_cyst] <- "Post-meiotic \ncyst"
seurat_integrated_ss$treatment <- "SR"
seurat_integrated_ss$treatment[grep("st", seurat_integrated_ss$sample)] <- "ST"
Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype

DefaultAssay(seurat_integrated_ss) <- "integrated"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

Xgenes <- filter(ortholog_table, chr == "Chr_X") %>% 
  dplyr::select(REF_GENE_NAME) %>% unlist() %>% 
  gsub("_", "-", .) %>% 
  c() %>% intersect(rownames(seurat_integrated_ss@assays$RNA@counts))
Agenes <- filter(ortholog_table, chr != "Chr_X") %>% 
  dplyr::select(REF_GENE_NAME) %>% unlist() %>% 
  gsub("_", "-", .) %>% 
  c() %>% intersect(rownames(seurat_integrated_ss@assays$RNA@counts))

DefaultAssay(seurat_integrated_ss) <- "RNA"


##################################################


sce <- seurat_integrated_ss %>% as.SingleCellExperiment(assay = "RNA")

### PSEUDOBULK INDIVIDUAL CELL TYPES ----
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]
y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))
yc <- calcNormFactors(y)

cpm <- edgeR::cpm(y, log=T) %>% as.data.frame()
cpmc <- edgeR::cpm(y, log=T) %>% as.data.frame()

heatmap(cor(cpm))
heatmap(cor(cpmc))

keep <- filterByExpr(y, group=agg_cell_types$treatment)
y <- y[keep,]

cpm <- edgeR::cpm(y, log=T) %>% as.data.frame()

raw_counts <- y$counts %>% as.data.frame()
cpm$genes <- rownames(cpm)
raw_counts$genes <- rownames(raw_counts)

cols <- c("treatment", "celltype", "sample")
tidy_cpm <- cpm %>% as.data.frame %>% 
  pivot_longer(cols = colnames(cpm)[1:(length(colnames(cpm))-1)], names_to = "sample", values_to = "logcpm") %>% 
  merge(y$samples[,cols], by.x = "sample", by.y = 0) %>% 
  rename(edgeRsample = sample,  sample = sample.y) %>% 
  mutate(Chr = ifelse(genes %in% Xgenes, "X", "A")) %>% 
  mutate(Chr = ifelse(!(genes %in% c(Xgenes, Agenes)), "Mito", Chr))

tidy_raw_counts <- raw_counts %>% as.data.frame %>% 
  pivot_longer(cols = colnames(raw_counts)[1:(length(colnames(raw_counts))-1)], names_to = "sample", values_to = "raw_counts") %>% 
  merge(y$samples[,cols], by.x = "sample", by.y = 0) %>% 
  rename(edgeRsample = sample,  sample = sample.y) %>% 
  mutate(Chr = ifelse(genes %in% Xgenes, "X", "A")) %>% 
  mutate(Chr = ifelse(!(genes %in% c(Xgenes, Agenes)), "Mito", Chr))

#################################################

check_dc_plot <- tidy_cpm %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic \ncyst", 
                                            "Post-meiotic \ncyst", "GSC & Spermatogonia", 
                                            "Spermatocytes", "Spermatids"))) %>%
  filter(Chr != 'Mito') %>% 
  group_by(celltype, genes, Chr) %>% 
  summarise(logcpm = mean(logcpm)) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = logcpm, fill = Chromosome)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "Cell Type", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 20) + 
  theme_classic() + scale_fill_brewer(palette = 'Set2')

check_dc_plot

ggsave("plots/dosage_compensation.pdf", check_dc_plot, width = 10, height = 6)

####################################################
