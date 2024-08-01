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

nFeatures_RNA_autosomes <- colSums(seurat_final@assays$RNA$counts[Agenes,] > 0)
seurat_final <- AddMetaData(seurat_final, nFeatures_RNA_autosomes, 'nFeature_RNA_autosomes')
X_exp <- PercentageFeatureSet(seurat_final, features = Xgenes)
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')

A_exp <- PercentageFeatureSet(seurat_final, features = Agenes)
seurat_final <- AddMetaData(seurat_final, A_exp, 'A_exp')

Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)/
  colSums(seurat_final@assays$RNA$counts > 1)
nXgene <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)
nXprop_toX <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)/
  colSums(seurat_final@assays$RNA$counts[Xgenes,])

seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')

seurat_final <- AddMetaData(seurat_final, X_exp < 10, 'Y_cells')
seurat_final <- AddMetaData(seurat_final, nXgene, 'nXgene')
seurat_final <- AddMetaData(seurat_final, nXprop_toX, 'nXprop_toX')

gene_names <- names(main_figure_markers)
mfms <- stringr::str_split(main_figure_markers, "\\(", simplify = T)[,1]
mfms[13] <- 'cup'
names(mfms) <- gene_names

gene_names_all_markers <- names(all_figure_markers)
afms <- stringr::str_split(all_figure_markers, "\\(", simplify = T)[,1]
afms[23] <- 'cup'
names(afms) <- gene_names_all_markers


##################################################

rm_axis = theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())


############## FIGURE 1 --------------------------
levels(seurat_final) <- c("Muscle", "Early cyst", "Late cyst", 
                          "GSC/Spermatogonia", "Primary spermatocytes", 
                          "Secondary spermatocytes", "Spermatids")
dps_all_figure <- DotPlot(seurat_final, features = names(mfms), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((mfms))) + 
  labs(y = '', x = "Marker gene expression") + 
  theme_classic() + 
  scale_colour_gradientn(colours = colorspace::diverge_hcl(8), 
                         labels = ) + 
  theme(legend.position="right", 
        legend.box = "horizontol") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title=element_text(size=13))



cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7")

nfeatures_figure <- seurat_final@meta.data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA_autosomes, fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + scale_fill_manual(values= cbPalette) + labs(
    x = "", y = "N\u00b0 detected autosomal genes") + 
  theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title=element_text(size=13))



# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

UMAP <- DimPlot(seurat_final, cols = cbPalette, group.by = 'celltype', label = T) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(legend.position="none")
UMAP2 <- UMAP[[1]]$data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  ggplot(aes(x = umap_1, y = umap_2, colour = celltype)) + 
  geom_point(stroke = NA, size = .7) + scale_color_manual(values= cbPalette) + 
  guides(colour = guide_legend(override.aes = list(size=3), title = "Cell type")) + 
  theme_classic() + labs(x = "UMAP 1", y = "UMAP 2", legend = 'Cell type') +
  theme(legend.position = 'none', 
        axis.title=element_text(size=15)) + 
  shadowtext::geom_shadowtext(data = UMAP[[1]]$layers[[2]]$data, aes(label = celltype), colour = 'black', 
                              size = 4, nudge_y = -0.5, nudge_x = 0, 
                              bg.color = 'white', bg.r = 0.1, bg.size = 0.5, 
                              fontface = "bold") + 
  geom_point(data = UMAP[[1]]$layers[[2]]$data, aes(x = umap_1, y = umap_2+0.1), colour = 'black',
             stroke = 1, size = 2, shape = 21, fill = 'white')



Fig1 <- 
  plot_grid(UMAP2, plot_grid(dps_all_figure + rm_axis, nfeatures_figure, nrow = 2, 
                                   align = 'v', labels = c("B", "C"), 
                                   label_x = -0.01), 
                  nrow = 1, align = 'h', rel_widths = c(1, 0.7, 0.7), axis = 't', 
                  labels =c("A", "")) 


ggsave("plots/FIG1_MarkersUMAPFeatures.pdf", Fig1, height = 9, width =14)
system("open plots/FIG1_MarkersUMAPFeatures.pdf")



S1a <- DimPlot(seurat_final, group.by = 'Phase') + 
  labs(title = "Cell-cycle Phase") + 
  labs(x = "UMAP 1", y = "UMAP 2")

S1b <- DotPlot(seurat_final, features = names(afms), assay = "RNA")+
  scale_x_discrete(labels = as.vector((afms))) + coord_flip() +
  labs(y = '', x = "Marker genes expression") + 
  theme_classic() + 
  scale_colour_gradientn(colours = colorspace::diverge_hcl(8), 
                         labels = ) + 
  theme(legend.position="right", 
        legend.box = "horizontol") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title=element_text(size=13))

S1c <- seurat_final@meta.data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA, fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + scale_fill_manual(values= cbPalette) + labs(
    x = "", y = "N\u00b0 detected genes") + 
  theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title=element_text(size=13))
  

FigS1 <- plot_grid(S1a, plot_grid(S1b + rm_axis, S1c, nrow = 2, 
                                   align = 'v', labels = c("B", "C"), 
                                   rel_heights = c(2.2,2)), 
                  nrow = 1, align = 'h', rel_widths = c(1, 0.7), axis = 't', 
                  labels =c("A", "")) 
ggsave("plots/FIGS1_MarkersUMAPFeatures.pdf", FigS1, height = 9, width =20)
system("open plots/FIGS1_MarkersUMAPFeatures.pdf")



############## FIGURE 2--------------------------
load("data/RData/DEG_DC.RData")

#test of normality 
check_normalities <- 
  XA %>% mutate(ratio = XA) %>% 
  filter(Chr == "X") %>% 
  filter(logcpm >= 1) %>% 
  group_by(celltype, treatment) %>% 
  summarise(p_value = shapiro.test(ratio)$p.value)

#Ratios not normally distributed 

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
  filter(Chr == "X") %>% 
  filter(logcpm >= 1) %>% 
  filter(treatment == "ST") %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = (XA), fill = treatment)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black') + 
  
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  rm_axis + theme(legend.position="right") + 
  guides(fill = 'none') +
  ylab(expression(log[2]~X:A~ratio)) + 
  geom_text(data = ST_wilcox_results, aes(x = celltype, y = 9, label = significance),
            vjust = -0.5, size = 5) + 
  ylim(-5,10)

expected_X_exp <- filter(ortholog_table, chr == "Chr_X") %>%
  dplyr::select(REF_GENE_NAME) %>% 
  unique() %>% nrow()
expected_A_exp <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2")) %>%
  dplyr::select(REF_GENE_NAME) %>% 
  unique() %>% nrow()
line = expected_X_exp/(expected_A_exp + expected_X_exp)


XCI_wilcox_results <- seurat_final@meta.data %>%
  filter(treatment == "ST") %>% 
  group_by(celltype, treatment) %>%
  summarise(
    p_value = wilcox.test(Xprop, mu = line, alternative = 'two.sided')$p.value
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.00001 ~ "***",
      p_value < 0.001 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ " "
    )
  )


XCI_plot <- seurat_final@meta.data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  filter(treatment == "ST") %>% 
  ggplot(aes(x = celltype, y = Xprop, fill = treatment)) + geom_violin() +
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  #geom_hline(yintercept = line, linetype = 'dashed', colour = 'black') +
  guides(fill = 'none') +
  #geom_text(data = XCI_wilcox_results, aes(x = celltype, y = 0.22, label = significance),
  #          vjust = -0.5, size = 5) +
  #ylim(0.05,0.25) +
  labs(y = "N\u00b0 expressed X genes / \n N\u00b0 expressed genes", x = "")



Genomic_activity <- cowplot::plot_grid(DC_plot, XCI_plot, ncol = 1, align="v", rel_heights = 
                                         c(1,1.2),labels = c("A", "B"),
                                       hjust = -0.5)
ggsave("plots/FIG2_Genomic_activity.pdf", Genomic_activity, height = 7, width = 4)
system("open plots/FIG2_Genomic_activity.pdf")



############## FIGURE 3 --------------------------
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
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + 
  stat_compare_means( aes(label = ..p.signif..), 
  label.x = 1.5, label.y = 9, method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
  symnum.args = list(
    cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
    symbols = c("***", "**", "*", " "))) + 
  
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black') + 
  
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  theme(legend.position="right") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab(expression(log[2]~X:A~ratio))
  

ggsave("plots/FIG3_SRvsST_DC.pdf", DC_plot_SR, height = 5, width = 5)
ggsave("plots/FIG3_SRvsST_DC.png", DC_plot_SR, height = 5, width = 5)

system("open plots/FIG3_SRvsST_DC.pdf")

#########################################


############## FIGURE 4 --------------------------

load("data/RData/tau_data.RData")

tau_plot_data <- tau_data %>% 
  mutate(tau_type = paste0(tau_type, " specific genes")) %>% 
  mutate(tau_type = factor(tau_type, levels = c("Muscle specific genes", "Cyst specific genes", 
                                                "Pre-meiosis specific genes", "Post-meiosis specific genes"))) %>% 
  mutate(chr = gsub("Chr_", "Chromosome  ", chr))

p1 <- tau_plot_data %>% filter(tau_type == "Pre-meiosis specific genes") %>%
  ggplot(aes(x = celltype, y = tau_exp, fill = treatment)) + geom_boxplot(outlier.alpha = 0.1) +
  labs(y = "% transcripts from \npre-meiosis-specific genes", x = "") + 
  scale_fill_brewer(palette = 'Set2') +
  facet_wrap(~chr, nrow = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme_classic() + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 0.58, method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
                      symnum.args = list(
                        cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
                        symbols = c("***", "**", "*", " "))) +
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  ylim(0,0.6) + rm_axis + 
  theme(strip.text.x = element_text(color = "white"),
        strip.background.x = element_rect(fill = "#4B4B4B", linetype = "solid"), 
        legend.position = c(0.92, 0.78))

p2 <- tau_plot_data %>% filter(tau_type == "Post-meiosis specific genes") %>%
  ggplot(aes(x = celltype, y = tau_exp, fill = treatment)) + geom_boxplot(outlier.alpha = 0.1) +
  labs(y = "% transcripts from \npost-meiosis-specific genes", x = "") + 
  scale_fill_brewer(palette = 'Set2') +
  facet_wrap(~chr, nrow = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme_classic() + 
  theme(legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  ylim(0,0.6) + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 0.58, method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
                      symnum.args = list(
                        cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
                        symbols = c("***", "**", "*", " "))) + 
  theme(strip.text.x = element_blank(),
        strip.background.x = element_rect(fill = "#4B4B4B", linetype = "solid"), 
        legend.position = 'none')
                                          

Fig4 <- ggarrange(plotlist = list(p1,p2), nrow = 2, common.legend = FALSE, heights = c(1,1.3), labels = c("(a)", "(b)"))
ggsave("plots/FIG4_tau_exp.pdf", Fig4, height = 8, width = 8)
system("open plots/FIG4_tau_exp.pdf")






