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
load("data/RData/ST_seurat_final.RData")
load("data/RData/ST_DEG_DC.RData")

DefaultAssay(seurat_final) <- "RNA"
DefaultAssay(ST_seurat_final) <- "RNA"

Idents(seurat_final) <- "celltype"
Idents(ST_seurat_final) <- "celltype"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

# Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
#   gsub("_", "-", .) %>% 
#   intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()


X_exp <- PercentageFeatureSet(seurat_final, features = Xgenes)
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')

A_exp <- PercentageFeatureSet(seurat_final, features = Agenes)
seurat_final <- AddMetaData(seurat_final, A_exp, 'A_exp')

Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)/
  colSums(seurat_final@assays$RNA$counts > 0)
# ST_Xgene_prop <- colSums(ST_seurat_final@assays$RNA$counts[Xgenes,] > 0)/
#   colSums(ST_seurat_final@assays$RNA$counts > 0)
nXgene <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)
nXprop_toX <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 0)/
  colSums(seurat_final@assays$RNA$counts[Xgenes,])

seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')
ST_seurat_final <- AddMetaData(object = ST_seurat_final, metadata = ST_Xgene_prop, col.name = 'Xprop')

seurat_final <- AddMetaData(seurat_final, X_exp < 10, 'Y_cells')
seurat_final <- AddMetaData(seurat_final, nXgene, 'nXgene')
seurat_final <- AddMetaData(seurat_final, nXprop_toX, 'nXprop_toX')

gene_names <- names(main_figure_markers)
mfms <- stringr::str_split(main_figure_markers, "\\(", simplify = T)[,1]
mfms[11] <- 'cup'
names(mfms) <- gene_names

##################################################

rm_axis = theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())


############## FIGURE 1 --------------------------
#seurat_final <- subset(ST_seurat_final, (treatment == "ST"))
seurat_final <- NormalizeData(seurat_final, normalization.method = "LogNormalize", scale.factor = 100)

levels(seurat_final) <- c("Muscle", "Pre-meiotic cyst", "Post-meiotic cyst", 
                          "GSC/Spermatogonia", "Primary Spermatocytes", 
                          "Secondary Spermatocytes", "Spermatids")
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

UMAP <- DimPlot(seurat_final, cols = cbPalette, group.by = 'celltype', label = T) + 
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
  theme(legend.position = 'none') + #c(0.97, 0.85)) + 
  shadowtext::geom_shadowtext(data = UMAP[[1]]$layers[[2]]$data, aes(label = celltype), colour = 'black', 
                              size = 4, nudge_y = -0.5, nudge_x = 0, 
                              bg.color = 'white', bg.r = 0.1, bg.size = 0.5, 
                              fontface = "bold") + 
  geom_point(data = UMAP[[1]]$layers[[2]]$data, aes(x = umap_1, y = umap_2), colour = 'black',
             stroke = 1, size = 2, shape = 21, fill = 'white')


#Arrange UMAP, Dotplot and nfeatures_figure with UMAP as first row and twice the height of the other two plots 
grid1 <- plot_grid(dps_all_figure, nfeatures_figure, nrow = 2, align = 'v', rel_heights = c(1,1), 
                   labels = c("B", "C"), hjust = 0.3)
Fig1 <- plot_grid(UMAP2, grid1, ncol = 2, align = 'v', rel_widths = c(1.4,1), 
                  labels = c("A"))
ggarrange(UMAP2, dps_all_figure, nfeatures_figure, nrow = 3)
plot_grid(UMAP2, nfeatures_figure, dps_all_figure, nrow = 1, align = 'h', 
          rel_widths = c(1, 1, 1.5), labels = c("A", "B", "C"))

ggsave("plots/FIG1_MarkersUMAPFeatures.pdf", Fig1, height = 8, width =14)


############## FIGURE 2--------------------------
load("data/RData/DEG_DC.RData")

XA <- ST_XA
#test of normality 
check_normalities <- 
  XA %>% mutate(ratio = XA) %>% 
  filter(Chr == "X") %>% 
  filter(logcpm > 1) %>% 
  filter(treatment == "ST") %>% 
  group_by(celltype) %>% 
  summarise(p_value = shapiro.test(ratio)$p.value)

#Ratios not normally distributed 

# Perform t-test for each group and store results in a data frame
ST_wilcox_results <- XA %>%
  filter(Chr == "X") %>% 
  filter(logcpm > 1) %>% 
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
      TRUE ~ "ns"
    )
  ) %>% 
  mutate(treatment = "ST")

DC_plot <- XA %>% 
  filter(Chr == "X") %>% 
  filter(logcpm > 1) %>% 
  filter(treatment == "ST") %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = (XA), fill = treatment)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black') + 
  
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  rm_axis + theme(legend.position="right") + 
  guides(fill = 'none') +
  ylab(expression(log[2](X:A)~ratio)) + 
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


XCI_wilcox_results <- ST_seurat_final@meta.data %>%
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
      TRUE ~ "ns"
    )
  )


XCI_plot <- ST_seurat_final@meta.data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  filter(treatment == "ST") %>% 
  ggplot(aes(x = celltype, y = Xprop, fill = treatment)) + geom_violin() +
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(yintercept = line, linetype = 'dashed', colour = 'black') +
  guides(fill = 'none') +
  geom_text(data = XCI_wilcox_results, aes(x = celltype, y = 0.22, label = significance),
            vjust = -0.5, size = 5) +
  ylim(0.05,0.25) +
  labs(y = "Proportion of expressed \ngenes X-linked", x = "")


XCI_plot2 <- ST_tidy_cpm_edger %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  group_by(celltype, sample, treatment) %>% 
  filter(Chr != "Mito") %>% 
  filter(treatment == "ST") %>% 
  summarise(Aexp = sum(logcpm > 0 & Chr == "A"), 
            Xexp = sum(logcpm > 0 & Chr == "X"),
            prop = Xexp/(Aexp + Xexp)) %>% 
  ggplot(aes(x = celltype, y = prop, fill = treatment)) + geom_boxplot() + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y = "Proportion of expressed \ngenes X-linked", x = "") + 
  geom_hline(yintercept = line, linetype = 'dashed', colour = 'black') +
  guides(fill = 'none')


Genomic_activity <- cowplot::plot_grid(DC_plot, XCI_plot, ncol = 1, align="v", rel_heights = 
                                         c(1,1),labels = c("A", "B"),
                                       hjust = -0.5)
ggarrange(DC_plot, XCI_plot, nrow = 2, labels = c("A", "B"))
ggsave("plots/FIG2_Genomic_activity.pdf", Genomic_activity, height = 7, width = 4)



############## FIGURE 3 --------------------------
DC_plot_SR <-  XA %>% 
  filter(Chr == "X") %>% 
  filter(logcpm > 1) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic cyst", 
                                                "Post-meiotic cyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = (XA), fill = treatment)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 12, method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
                      symnum.args = list(
                        cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
                        symbols = c("***", "**", "*", "ns"))) + 
  
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black') + 
  
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  theme(legend.position="right") + 
  #guides(fill = 'none') +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab(expression(log[2](X:A)~ratio))


ggsave("plots/FIG3_SRvsST_DC.pdf", DC_plot_SR, height = 5, width = 5)


#########################################


############## FIGURE 4 --------------------------



