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

#### LOAD DATA AND PREP DATA ---
load("data/RData/seurat_final.RData")
DefaultAssay(seurat_final) <- "RNA"
Idents(seurat_final) <- "celltype"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()


X_exp <- PercentageFeatureSet(seurat_final, features = Xgenes)
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')

A_exp <- PercentageFeatureSet(seurat_final, features = Agenes)
seurat_final <- AddMetaData(seurat_final, A_exp, 'A_exp')

Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)/
  colSums(seurat_final@assays$RNA$counts > 0)
nXgene <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)
nXprop_toX <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)/
  colSums(seurat_final@assays$RNA$counts[Xgenes,])

seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')
seurat_final <- AddMetaData(seurat_final, X_exp < 10, 'Y_cells')
seurat_final <- AddMetaData(seurat_final, nXgene, 'nXgene')
seurat_final <- AddMetaData(seurat_final, nXprop_toX, 'nXprop_toX')

gene_names <- names(main_figure_markers)
mfms <- stringr::str_split(main_figure_markers, "\\(", simplify = T)[,1]
mfms[11] <- 'cup'
names(mfms) <- gene_names

##################################################


############## FIGURE 1 --------------------------
dps_all_figure <- DotPlot(seurat_final, features = names(mfms), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((mfms))) + 
  labs(y = '', x = "Genes") + 
  theme_classic() + 
  scale_colour_gradientn(colours = colorspace::diverge_hcl(8), 
                         labels = ) + 
  #theme(legend.position="right") + 
  theme(legend.position="right") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
  
  



cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7")

nfeatures_figure <- seurat_final@meta.data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA, fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + scale_fill_manual(values= cbPalette) + labs(
    x = "", y = "Number of expressed genes") + 
  theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))



# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

UMAP <- DimPlot(seurat_final, cols = cbPalette, group.by = 'celltype') + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(legend.position="none")
UMAP2 <- UMAP[[1]]$data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  ggplot(aes(x = umap_1, y = umap_2, colour = celltype)) + 
  geom_point(stroke = NA, size = .7) + scale_color_manual(values= cbPalette) + 
  guides(colour = guide_legend(override.aes = list(size=3), title = "Cell type")) + 
  theme_classic() + 
  labs(x = "UMAP 1", y = "UMAP 2", legend = 'Cell type') +
  theme(legend.position = c(0.97, 0.85)) + 
  xlim(c(-15, 15)) + ylim(c(-15, 15))
UMAP2
#Arrange UMAP, Dotplot and nfeatures_figure with UMAP as first row and twice the height of the other two plots 
grid1 <- plot_grid(dps_all_figure, nfeatures_figure, nrow = 2, align = 'v', rel_heights = c(1,1), 
                   labels = c("B", "C"), hjust = 0.3)
Fig1 <- plot_grid(UMAP2, grid1, ncol = 2, align = 'v', rel_widths = c(1.4,1), 
          labels = c("A"))

ggsave("plots/FIG1_MarkersUMAPFeatures.pdf", Fig1, height = 8, width =14)


############## FIGURE 2--------------------------

tidy_cpm_all <- read.csv("data/MANUSCRIPT/all_cpm_values.csv")
tidy_cpm_filt <- read.csv("data/MANUSCRIPT/filtered_cpm_values.csv")

rm_axis = theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())

DC_plot <- tidy_cpm_filt %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  filter(treatment == "ST" & Chr != 'Mito' & logcpm > 1) %>% 
  group_by(celltype, genes, Chr, treatment) %>% 
  summarise(logcpm = log2(mean(2^logcpm))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = logcpm, fill = Chromosome)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  stat_compare_means( aes(label = ..p.signif..), 
  label.x = 1.5, label.y = 20, method = 't.test') + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  rm_axis + theme(legend.position="right")

expected_X_exp <- filter(ortholog_table, chr == "Chr_X") %>%
  dplyr::select(REF_GENE_NAME) %>% 
  unique() %>% nrow()
expected_A_exp <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2")) %>%
  dplyr::select(REF_GENE_NAME) %>% 
  unique() %>% nrow()
line = expected_X_exp/(expected_A_exp + expected_X_exp)

XCI_plot <- tidy_cpm_filt %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                          "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                          "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                          "Spermatids"))) %>% 
  group_by(celltype, sample, treatment) %>% 
  filter(Chr != "Mito") %>% 
  filter(treatment == "ST") %>% 
  summarise(Aexp = sum(logcpm > 1 & Chr == "A"), 
            Xexp = sum(logcpm > 1 & Chr == "X"),
            prop = Xexp/(Aexp + Xexp)) %>% 
  ggplot(aes(x = celltype, y = prop)) + geom_boxplot() + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y = "Proportion of expressed \ngenes X-linked", x = "") + 
  geom_hline(yintercept = line, linetype = 'dashed', colour = 'black')

XCI_plot

Genomic_activity <- cowplot::plot_grid(DC_plot, XCI_plot, ncol = 1, align="v", rel_heights = 
                                         c(1,1),labels = c("A", "B"),
                                       hjust = -0.5)

ggsave("plots/FIG2_Genomic_activity.pdf", Genomic_activity, height = 7, width = 6)




################ DC in SR Plot --------------
DC_plot_SR <- tidy_cpm_filt %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  filter(treatment == "SR") %>%
  filter(Chr != 'Mito') %>% 
  filter(logcpm > 1) %>% 
  group_by(celltype, genes, Chr, treatment) %>% 
  summarise(logcpm = log2(mean(2^logcpm))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = logcpm, fill = Chromosome)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  stat_compare_means( aes(label = ..p.signif..), 
                     label.x = 1.5, label.y = 20, method = 't.test') + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
   theme(legend.position="none")

STSR_DC <- cowplot::plot_grid(DC_plot, DC_plot_SR, ncol = 1, align="v", rel_heights = 
                                         c(1,1.5),labels = c("A", "B"),
                                       hjust = -0.5)
ggsave("plots/FIGSX_SRvsST_DC.pdf", STSR_DC, height = 7, width = 6)




######## delete ----
tidy_cpm2 %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  group_by(celltype, sample, treatment) %>% 
  filter(Chr != "Mito") %>% 
  filter(treatment == "ST") %>% 
  summarise(Aexp = sum(logcpm > 2 & Chr == "A"), 
            Xexp = sum(logcpm > 2 & Chr == "X"),
            prop = Xexp/(Aexp + Xexp)) %>% 
  ggplot(aes(x = celltype, y = prop, fill = treatment)) + geom_boxplot() + theme_classic()


sr1_cols <- which(seurat_final@meta.data$celltype == "Muscle" & 
                    seurat_final@meta.data$sample == "sr1")

seurat_final@assays$RNA@counts[,sr1_cols] %>% 
  rowSums() %>% .['PB.1']

Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)/
  colSums(seurat_final@assays$RNA$counts[c(Xgenes, Agenes),] > 0)

Xexp <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)
Aexp <- colSums(seurat_final@assays$RNA$counts[Agenes,] > 0)
xgp_df <- data.frame(Xprop = Xgene_prop, bc = names(Xgene_prop), Xexp = Xexp, Aexp = Aexp)
xgp_df <- merge(xgp_df, seurat_final@meta.data[,c("sample", "treatment", "celltype", "cells")], 
                by.x = 'bc', by.y = 'cells')

xgp_df %>% group_by(celltype, sample, treatment) %>% 
  summarise(prop = mean(Xprop)) %>% 
  ggplot(aes(x = celltype, y = prop, fill = treatment)) + 
  geom_boxplot() 

seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = Xprop)) + geom_boxplot() + 
  facet_wrap(~sample)






