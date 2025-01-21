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
load("data/trajectory/sce_GAMed_8k.RData")
load("data/RData//traj_analyis.RData")

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

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

## MARKER GENES ##
markers <- read.xlsx("data/MANUSCRIPT/Drive Manuscript Supplementary Tables.xlsx", sheet = 2)[1:30,c(1:6)] 
row_order <- c("M", "CYSC", "EC", "LC", "C", "E", "G", "G, PS", "PS, SS", "PS, SS, ST", "PS, ST", "ST", "LST")
markers$Drosophila.cell.type <- factor(markers$Drosophila.cell.type, levels = row_order)
markers <- markers[order(markers$Drosophila.cell.type),]

mk_genes <- markers$Teleopsis.dalmanni.Ortholog %>% gsub(",", " ", .) %>% 
  str_split(" ") %>% unlist() %>% 
  gsub("_", "-", .) %>% intersect(rownames(seurat_final@assays$RNA))


new_name_func <- function(gene){
  nn <- paste(markers$`Drosophila.Marker.(Gene.Name)`[
    grep(gene, gsub("_", "-", markers$Teleopsis.dalmanni.Ortholog ))], 
    collapse = ", ")
  nn2 <- paste0(nn , " (", gene, ")")
  return(nn)
}

new_names <- sapply(mk_genes, new_name_func)
#duplicate genes
rm <- c("PB.4544", "PB.4546", "gene−9679", "STRG.14307", "g9515")
all_figure_markers <- new_names[!names(new_names) %in% rm]
all_figure_markers <- all_figure_markers[-which(duplicated(all_figure_markers))]

main_figure_markers <- all_figure_markers[c(1,3:7,14,15,18,21,23,25,24,26,28)]


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


### TITLE CHANGES FOR PLOTS ###
gene_names <- names(main_figure_markers)
mfms <- stringr::str_split(main_figure_markers, "\\(", simplify = T)[,1]
mfms[14] <- 'cup'
names(mfms) <- gene_names

gene_names_all_markers <- names(all_figure_markers)
afms <- stringr::str_split(all_figure_markers, "\\(", simplify = T)[,1]
afms[26] <- 'cup'
names(afms) <- gene_names_all_markers

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

dps_all_figure <- DotPlot(seurat_final, features = names(mfms), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((mfms))) + 
  labs(y = '', x = "Marker gene expression") + 
  theme_classic() + scale_colour_gradientn(colours = colorspace::diverge_hcl(8)) + 
  theme(legend.position="right", legend.box = "vertical", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) + 
  guides(color = guide_colorbar(title = 'Scaled expression'), 
        size = guide_legend(title = '% of cells expressed'))
  


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
                              size = 4, nudge_y = -0.5, nudge_x = 0, 
                              bg.color = 'white', bg.r = 0.1, bg.size = 0.5, 
                              fontface = "bold") + 
  geom_point(data = UMAP[[1]]$layers[[2]]$data, aes(x = umap_1, y = umap_2+0.1), colour = 'black',
             stroke = 1, size = 2, shape = 21, fill = 'white')

# DC Plot #
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
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
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

traj_umap <- ggplot() + 
  geom_point(data = Pseudotime_data, aes(x = umap_1, y = umap_2, colour = Pseudotime), alpha = 0.3, stroke = NA) + 
  geom_path(data = traj_line[[1]], aes(x = umap_1, y = umap_2), colour = 'black', size = 0.8, alpha = 0.7) +
  theme_classic() +
  scale_color_continuous(low = "purple", high = "orange") + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(legend.position = c(0.8, 0.75))

traj_nfeatures <- Pseudotime_data[Pseudotime_data$treatment == "ST",]  %>% 
  ggplot(aes(x = Pseudotime, y = nFeature_RNA_autosomes, colour = `Cell type`)) + 
  geom_point(alpha = 0.2, stroke=NA) + scale_colour_manual(values= cbPalette) +
  theme_classic() + 
  geom_line(stat = 'smooth', aes(x =Pseudotime, y = nFeature_RNA_autosomes, colour = NA),  colour = 'black', method = 'loess', se = F, 
            alpha = 0.7, size = 1) + 
  labs(y = "N\u00b0 expressed \nautosomal genes", x = "") + 
  theme(legend.position = 'none', 
        legend.background=element_blank()) + rm_axis


traj_XoverAll <- Pseudotime_data[Pseudotime_data$treatment == "ST",] %>% 
  ggplot(., aes(x = Pseudotime, y = Xprop, colour = `Cell type`)) + 
  geom_point(alpha = 0.2, stroke = NA) + scale_colour_manual(values= cbPalette) +
  geom_line(stat = 'smooth', aes(x =Pseudotime, y = Xprop, colour = NA),  colour = 'black', method = 'loess', se = F, 
            alpha = 0.7, size = 1) + 
  theme_classic() + 
  theme(legend.position = 'none', 
        legend.background=element_blank()) +
  labs(y ="N\u00b0 expressed X genes / \n N\u00b0 total expressed genes", x = "") + rm_axis

traj_XoverAll_SR <- Pseudotime_data %>% 
  mutate(Treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  ggplot(., aes(x = Pseudotime, y = Xprop)) + 
  #geom_point(alpha = 0.01, aes(colour = `Cell type`)) + scale_colour_manual(values= cbPalette) +
  guides(colour = 'none') + 
  ggnewscale::new_scale_color() + 
  geom_smooth(method = 'loess', alpha = 0.5, size = 1, 
              aes(color = Treatment, fill = Treatment, ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  guides(colour= guide_legend(override.aes=list(fill=NA))) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = 'Set2') +
 theme_classic() + 
 theme(legend.position = 'top', 
      legend.background=element_blank()) +
  
 labs(y ="N\u00b0 expressed X genes / \n N\u00b0 total expressed genes", x = "Pseudotime", colour = "", fill = "") + 
  rm_axis

 
 traj_celltypes <- Pseudotime_data[Pseudotime_data$treatment == "ST",] %>% 
  ggplot(aes(x = Pseudotime, y = `Cell type`, fill = `Cell type`)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  theme_classic() + theme(legend.position = 'none') +
  scale_fill_manual(values= cbPalette) + 
  labs(y = "") + 
  scale_y_discrete(limits=rev) + theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
                                       axis.ticks = element_line(color = "black"))


DC_plot_SR <-  XA %>% 
  filter(Chr == "X") %>% 
  filter(logcpm >= 1) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = (XA), fill = treatment)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "") + 
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
  ylab(expression(log[2]~X:A~expression~ratio)) + 
  scale_y_continuous(breaks=c(-5,0,5)) 


dif_exp_figure <- dif_exp_data %>%
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
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


#Figure 1 
A <- ggplotGrob(UMAP_final)
B <- ggplotGrob(dps_all_figure)
Fig1 <- plot_grid(A,B, labels = c("(a)", "(b)"), align = 'v', axis = 'l', ncol = 1, rel_heights = c(8.5,7), label_x = -.015)
ggsave("plots/FIG1.pdf", Fig1, width =8, height = 11)
system("open plots/FIG1.pdf")

traj_umap <- ggplot() + 
  geom_point(data = Pseudotime_data, aes(x = umap_1, y = umap_2, colour = Pseudotime), alpha = 0.3, stroke = NA) + 
  geom_path(data = traj_line[[1]], aes(x = umap_1, y = umap_2), colour = 'black', size = 0.8, alpha = 0.7) +
  theme_classic() +
  scale_color_continuous(low = "purple", high = "orange") + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(legend.position = c(0.85, 0.28)) + 
  ylim(-12,12)

traj_umap
#Figure 2
A <- (traj_nfeatures)
B <- (traj_celltypes)
C <- (traj_umap)
AB <- plot_grid(A, B, align = 'v', axis = 'l', nrow = 2, labels = c("(a)", "(b)"), rel_heights = c(3,2), label_size = 10.5)
Fig2 <- plot_grid(AB, C,  ncol = 2, labels = c("", "(c)"), rel_widths = c(2,2), label_size = 10.5)
ggsave("plots/FIG2.pdf", Fig2, width = 7, height = 3.5)
system("open plots/FIG2.pdf")


#Figure 3
A <- ggplotGrob(traj_XoverAll)
B <- ggplotGrob(traj_celltypes)
C <- ggplotGrob(DC_plot)
BC <- align_plots(B, C, align = 'h', axis = 'b')
AB <- plot_grid(A,BC[[1]], align = 'v', axis = 'l', nrow = 2, labels = c("(a)", "(b)"), rel_heights = c(3,3), label_size = 10.5)
Fig3 <- plot_grid(AB, BC[[2]], ncol = 2, labels = c("", "(c)"), rel_widths = c(2,2), label_size = 10.5)
ggsave("plots/FIG3.pdf", Fig3, width = 7, height = 4.5)
system("open plots/FIG3.pdf")


#Figure 4
A1 <- ggplotGrob(traj_XoverAll_SR)
A2 <- ggplotGrob(traj_celltypes)
B <- ggplotGrob(DC_plot_SR)
C <- ggplotGrob(dif_exp_figure)

A2BC <- align_plots(A2, B,C, align = "h", axis = "b")
A1BC <- align_plots(A1, A2BC[[2]], A2BC[[3]], align = "h", axis = "t")
A <- plot_grid(A1BC[[1]], A2BC[[1]], ncol = 1, rel_heights = c(4,4), align = 'v', axis = 'l', labels = c("", "(b)"))
Fig4 <- plot_grid(A, A2BC[[2]], A2BC[[3]], ncol = 3, labels = c("(a)", "(c)", "(d)"))
ggsave("plots/FIG4.pdf", Fig4, width = 10, height = 4.5)
system("open plots/FIG4.pdf")







##################  FIGURE S1 #######################

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7", "darkgrey")

S1a <- DimPlot(seurat_final, group.by = 'Phase') + 
  labs(title = "") + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(axis.text.x = element_text(color="black"), 
        axis.ticks = element_line(color = "black"), 
        legend.position = c(0.1, 0.8))

S1b <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA, fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + 
  scale_fill_manual(values= cbPalette) + 
  labs(x = "", y = "N\u00b0 expressed genes") + 
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black"))


S1c <- DotPlot(seurat_final, features = names(afms), assay = "RNA")+
  scale_x_discrete(labels = as.vector((afms))) + coord_flip() +
  labs(y = '', x = "Marker genes expression") + 
  theme_classic() + 
  scale_colour_gradientn(colours = colorspace::diverge_hcl(8)) + 
  theme(legend.position="right", 
        legend.box = "horizontol", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black"))

SA <- ggplotGrob(S1a)
SB <- ggplotGrob(S1b)
SC <- ggplotGrob(S1c)

##
AB <- plot_grid(SA, SB, align = 'v', axis = 'l', nrow = 2, labels = c("(a)", "(b)"), rel_heights = c(3,3))
FigS1 <- plot_grid(AB, SC, ncol = 2, labels = c("", "(c)"), rel_widths = c(2,3))
ggsave("plots/S1.pdf", FigS1, height = 7, width =10)
system("open plots/S1.pdf")

#####################################################

##################  FIGURE S2 #######################

#PLOIDY PLOT IN SCRIPT Ploidy_check.R

#####################################################


##################  FIGURE S3 #######################


cts <- seurat_final@meta.data$celltype %>% levels()
comparisons <- list(cts[4:5], cts[5:6], cts[6:7], cts[7:8])

XCI_plot <- seurat_final@meta.data %>% 
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
  ylim(0,0.29) + 
  rm_axis

No_x_expressed_plot <- seurat_final@meta.data %>% 
  mutate(treatment = factor(treatment, levels = c("ST", "SR"))) %>% 
  filter(treatment == "ST") %>% 
  ggplot(aes(x = celltype, y = Xprop2, fill = treatment)) + geom_boxplot(outlier.alpha = 0.1) +
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  theme(legend.position="none", legend.box = "vertical", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) +
  labs(y = "N\u00b0 expressed X genes / \n N\u00b0 X genes", x = "")

S3 <- plot_grid(XCI_plot, No_x_expressed_plot, align = 'v', axis = 'l', nrow = 2, labels = c("(a)", "(b)"), rel_heights = c(5,5))

ggsave("plots/FIGS3_XCI.pdf", S3, height = 7, width = 5)
system("open plots/FIGS3_XCI.pdf")

#####################################################


#####################################

#Find scripts for S4 in Data_Prep.R

#####################################




##################  FIGURE S5 #######################

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7", "darkgrey")

S5_data <- ortholog_table %>%
  group_by(REF_GENE_NAME) %>%
  summarise(start = min(start), end = max(end), REF_GENE_NAME = REF_GENE_NAME[1], 
            consensus_gene = consensus_gene[1], chromosome = chr[1]) %>%
  merge(XA, by.x = "REF_GENE_NAME", by.y = "genes") %>% 
  mutate(chromosome = gsub("Chr_", "Chromosome ", chromosome)) %>%
  mutate(celltype = gsub(" sp", "\n sp", celltype)) %>% 
  mutate(celltype = gsub("/", "/\n", celltype)) %>% 
  filter(treatment == 'ST') %>% 
   mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                 "Late cyst", "GSC/\nSpermatogonia", 
                                                 "Primary\n spermatocytes", "Secondary\n spermatocytes", 
                                                "Early\n spermatids", "Late\n spermatids")))

S5a_lines <- S5_data %>% 
  group_by(celltype, chromosome) %>% 
  summarise(cpm = log2(median(2^logcpm))) 

comps <- list(c("Chromosome 1", "Chromosome 2"), c("Chromosome 2", "Chromosome X"), c("Chromosome 1", "Chromosome X"))

a <- S5_data %>%
  ggplot(aes(y = logcpm, fill = chromosome, x = celltype)) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = 'Dark2') +
  theme_classic() +
  theme(legend.position = 'top', axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) + 
  labs(fill = "", x = "", y = log[2]~"(CPM)") + 
  geom_pwc(method = 'wilcox.test', hide.ns = T, label = "p.signif",step.increase = 0.07, tip.length = 0.005, 
           vjust = 0.5, symnum.args = list(
                                 cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
                                    symbols = c("***", "**", "*", "ns")))
  
b <- S5_data %>% 
  ggplot(aes(x = start/1000000, y = logcpm, colour = chromosome)) + geom_point(alpha = 0.1) + 
  theme_classic() +
  facet_grid(celltype ~ chromosome, switch = "both", scales = 'free_x') +
  theme(strip.placement = "outside") + 
  ylab(log[2]~"(CPM)") + xlab("Chromosomal position (Mb)") + 
  theme(strip.background = element_blank(), 
        legend.position = 'none') + 
  geom_hline(data = S5a_lines, aes(yintercept = cpm), colour = 'black') + 
  scale_color_brewer(palette = 'Dark2')
  
S5a <- ggplotGrob(a)
S5b <- ggplotGrob(b)
S5 <- plot_grid(S5a, S5b, ncol = 1, labels = c("(a)", "(b)"), rel_heights = c(2,5), align = 'v', axis = 'l')


ggsave("plots/FIGS5_chr_exp.pdf", S5, height = 12, width = 7)
system("open plots/FIGS5_chr_exp.pdf")


#####################################################


############## FIGURE S6 TAU DATA ###############
# Plotting tau data 
tau_data <- read.csv("data/MANUSCRIPT/tau_data.csv")

filtered_tau_data <- tau_data %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  filter(ts %in% c("Testes", "Universal", "Weak testes")) %>% 
  mutate(ts = factor(ts, levels = c("Testes", "Weak testes", "Universal"))) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids")))


tau_wilcox_results <- filtered_tau_data %>%
  group_by(celltype, ts) %>%
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
  )

tau_plots <- filtered_tau_data %>% 
  ggplot(aes(x = celltype, y = XA, fill = ts)) + geom_boxplot() + 
  labs(x = "", y = "log(CPM)", fill = "Tau class") + 
  
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black') + 
  theme_classic() + scale_fill_brewer(palette = 'Paired') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  theme(legend.position="top", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) +
  ylab(expression(log[2]~X:A~expression~ratio)) + 
  labs(fill = 'Tissue-specificity') +
  scale_y_continuous(breaks=c(-5,0,5)) + 
  geom_text(data = tau_wilcox_results[tau_wilcox_results$ts == "Testes",], aes(x = celltype, y = 6, label = significance),
            nudge_x = -0.18,vjust = -0.5, size = 5) + 
  geom_text(data = tau_wilcox_results[tau_wilcox_results$ts == "Universal",], aes(x = celltype, y = 7, label = significance),
            nudge_x = 0.18,vjust = -0.5, size = 5)



ggsave("plots/S6_tau_XA_plot.pdf", tau_plots, height = 5, width = 5)
system("open plots/S6_tau_XA_plot.pdf")

#####################################################




##################  FIGURE S7 #######################

#See MSL_plots.R

#####################################################




##################  FIGURE S8 #######################

XCI_plot <- seurat_final@meta.data %>% 
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

S8 <- ggplotGrob(XCI_plot)
ggsave("plots/FIGS8_XCI.pdf", S8, height = 4, width = 5)
system("open plots/FIGS8_XCI.pdf")


#####################################################



##################  FIGURE S9 #######################
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7", "darkgrey")

myTheme <- theme(legend.text = element_text(size = 10), 
                 legend.title = element_text(size = 12),
                 legend.key.size = unit(1, 'line'))

dif_exp_data <- dif_exp_data %>% 
  mutate(Chromosome = chr) %>% 
  mutate(`Expression bias` = Significant) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids")))


dif_exp_data$Chromosome[dif_exp_data$Chromosome == "Chr_X"] <- "X"
dif_exp_data$Chromosome[dif_exp_data$Chromosome == "Chr_1"] <- "Autosomes"
dif_exp_data$Chromosome[dif_exp_data$Chromosome == "Chr_2"] <- "Autosomes"
labels <- dif_exp_data %>% 
  filter(abs(logFC) > 2 & FDR < 0.05 & consensus_gene != genes) %>% 
  filter(!grepl("CG", consensus_gene)) %>%
  group_by(celltype) %>% 
  top_n(FDR, n = -12)

strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = cbPalette, 
                                                                  alpha = 0.5))
DIF_EXP_FIGa <- dif_exp_data %>% as.tibble() %>% 
  filter(Chromosome == "X") %>% 
  mutate(`Expression bias` = factor(`Expression bias`, levels = c("ST-biased", "SR-biased",  
                                                                  "Unbiased"))) %>% 
  ggplot(aes(x = logFC, y = -log10(FDR), 
             colour = `Expression bias`)) +
  geom_point() +
  scale_colour_brewer(palette = 'Set2') +
  geom_hline(yintercept = -log10(0.05), linetype = 2) + 
  geom_vline(xintercept = c(-1, 1), linetype = 2) + 
  theme_classic() +
  ggh4x::facet_wrap2(~ celltype, strip = strip, nrow = 2) + 
  theme(#strip.text.x = element_text(color = "white"),
    #strip.background.x = element_rect(fill = "#4B4B4B", linetype = "solid"), 
    legend.position = "none") + 
  
  geom_label_repel(data = labels[labels$Chromosome == 'X',], aes(x = logFC, y = -log10(FDR), label = consensus_gene),
                   box.padding   = 0.9, 
                   point.padding = 0.1, force = 120,
                   segment.color = 'grey50',show.legend = F, size = 2.8) +
  xlim(-7.5, 7.5) + 
  ylim(0,10) + myTheme + 
  guides(shape = 'none', 
         colour = guide_legend(override.aes = list(size = 2))) + 
  labs(x = expression(log[2]~fold~change), 
       y= expression(-log[10]~(FDR)))


DIF_EXP_FIGb <- dif_exp_data %>% as.tibble() %>% 
  filter(Chromosome != "X") %>% 
  mutate(`Expression bias` = factor(`Expression bias`, levels = c("ST-biased", "SR-biased",  
                                                                  "Unbiased"))) %>% 
  ggplot(aes(x = logFC, y = -log10(FDR), 
             colour = `Expression bias`)) +
  geom_point() +
  scale_colour_brewer(palette = 'Set2') +
  geom_hline(yintercept = -log10(0.05), linetype = 2) + 
  geom_vline(xintercept = c(-1, 1), linetype = 2) + 
  theme_classic() +
  ggh4x::facet_wrap2(~ celltype, strip = strip, nrow = 2) + 
  theme(#strip.text.x = element_text(color = "white"),
    #strip.background.x = element_rect(fill = "#4B4B4B", linetype = "solid"), 
    legend.position = "bottom") + 
  
  geom_label_repel(data = labels[labels$Chromosome != "X",], aes(x = logFC, y = -log10(FDR), label = consensus_gene),
                   box.padding   = 0.3, 
                   point.padding = 0.3, force =100,
                   segment.color = 'grey50',show.legend = F, size = 2.8) +
  xlim(-7.5, 7.5) + 
  ylim(0,10) + myTheme + 
  labs(x = expression(log[2]~fold~change), 
       y= expression(-log[10]~(FDR)), 
       colour = '')



S9 <- ggarrange(DIF_EXP_FIGa, DIF_EXP_FIGb, ncol = 1, nrow = 2, labels = c("(a)", "(b)"), heights = c(1,1.1))
ggsave("plots/FIGS9_dif_exp.pdf", S9, height = 11, width = 10)
system("open plots/FIGS9_dif_exp.pdf")

#####################################################



##################  FIGURE S10 #######################


S10_dif_exp_data <- dif_exp_data %>% 
  mutate(Significant = factor(Significant, levels = c("ST-biased", "SR-biased", "Unbiased"))) %>% 
  mutate(alpha = ifelse(Significant != "Unbiased", 1, 0.001)) %>% 
  mutate(chr = gsub("Chr_", "Chromosome ", chr)) %>%
  mutate(celltype = gsub(" sp", "\n sp", celltype)) %>% 
  mutate(celltype = gsub("/", "/\n", celltype)) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/\nSpermatogonia", 
                                                "Primary\n spermatocytes", "Secondary\n spermatocytes", 
                                                "Early\n spermatids", "Late\n spermatids"))) %>% 
  ggplot(aes(x = start/1000000, y = logFC, alpha = alpha,  colour = Significant)) + geom_point() +
  theme_classic() + 
  facet_grid(celltype ~ chr, switch = "both", scales = 'free_x') +
  theme(strip.placement = "outside") +
  ylab(expression(log[2]~fold~change)) + xlab("Chromosomal position (Mb)") + 
  scale_color_brewer(palette = 'Set2') + 
  theme(strip.background = element_blank()) + 
  guides(alpha = 'none') + labs(colour = "Expression bias")


ggsave("plots/FIGS10_chr_dif_exp.pdf", S10_dif_exp_data, height = 10, width = 8)
system("open plots/FIGS10_chr_dif_exp.pdf")

#####################################################


######## S11 #############
#See Trajectory_Analaysis.R for S11 plot
##########################


