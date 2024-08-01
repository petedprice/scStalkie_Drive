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
#rm(list = ls())

#### LOAD DATA AND PREP DATA ---
#load("data/RData/seurat_final.RData")

#DefaultAssay(seurat_final) <- "RNA"


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
  theme(legend.position="bottom", 
        legend.box = "vertical") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
    



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
  theme(legend.position="none") + rm_axis + 
  theme(axis.line.x = element_blank())
# theme(axis.text.x = element_text(angle = 45, hjust=1))#, 
#        axis.title=element_text(size=13))



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
  theme(legend.position = 'none') + #, 
#        axis.title=element_text(size=15)) + 
  shadowtext::geom_shadowtext(data = UMAP[[1]]$layers[[2]]$data, aes(label = celltype), colour = 'black', 
                              size = 4, nudge_y = -0.5, nudge_x = 0, 
                              bg.color = 'white', bg.r = 0.1, bg.size = 0.5, 
                              fontface = "bold") + 
  geom_point(data = UMAP[[1]]$layers[[2]]$data, aes(x = umap_1, y = umap_2+0.1), colour = 'black',
             stroke = 1, size = 2, shape = 21, fill = 'white')



DC_plot <- XA %>% 
  filter(Chr == "X") %>% 
  filter(logcpm >= 1) %>% 
  filter(treatment == "ST") %>% 
  #mutate(celltype = recode(celltype, "P\u00b0 spermatocytes" = "Primary spermatocytes", 
  #                         "S\u00b0 spermatocytes" = "Secondary spermatocytes")) %>%
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = (XA), fill = treatment)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black', linewidth = 0.7) + 
  
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  #rm_axis + theme(legend.position="right") + 
  guides(fill = 'none') +
  ylab(expression(log[2]~X:A~ratio)) + 
  geom_text(data = ST_wilcox_results, aes(x = celltype, y = 9, label = significance),
            vjust = -0.5, size = 5) + 
  ylim(-5,10) + rm_axis + 
  theme(axis.line.x = element_blank())

#  theme(axis.text.x = element_text(angle = 45, hjust=1))#, 
#        axis.title=element_text(size=13))

A <- ggplotGrob(UMAP2)
B <- ggplotGrob(dps_all_figure)
C <- ggplotGrob(nfeatures_figure)
D <- ggplotGrob(DC_plot)
E <- ggplotGrob(XCI_plot)
AD <- align_plots(A,D, align = 'h', axis = 'b')
CD <- plot_grid(C, AD[[2]], labels = c("C", "D"), align = 'v', ncol = 1, rel_heights = c(1, 1.5))
ACD <- plot_grid(AD[[1]], CD, labels = c("A", ""), rel_widths = c(1, 0.5))
Fig1 <- plot_grid(ACD, B, nrow = 2, axis = 'l', rel_heights = c(1, 0.5), labels = c("", "B"))
ggsave("plots/FIG1.pdf", Fig1, width = 13, height = 13)
#system("open plots/FIG1.pdf")

AC <- align_plots(A,C, align = 'h', axis = 'b')
BE <- align_plots(B,E, align = 'h', axis = 'b')
AB <- plot_grid(AC[[1]],BE[[1]], align = 'v', nrow = 2, rel_heights = c(3, 2), 
                labels = c("(a)", "(b)"))

CDE <- plot_grid(AC[[2]],D,BE[[2]], align = 'v', axis = 'b', nrow = 3, rel_heights = c(1,1,1.6), 
                 labels = c("(c)", "(d)", "(e)"))

Fig1 <- plot_grid(AB, CDE, align = 'h', ncol = 2, rel_widths = c(2, 1))
ggsave("plots/FIG1.pdf", Fig1, width = 13, height = 13)
system("open plots/FIG1.pdf")


###############################
AB <- align_plots(A,B, align = 'v', axis = 'l')
DE <- align_plots(D,E, align = 'v', axis = 'l')

BE <- plot_grid(AB[[2]],DE[[2]], align = 'h', axis = 'bt', labels = c("(b)", "(e)"), rel_widths = c(2,1))
CD <- plot_grid(C,DE[[1]], align = 'v', axis = 'l', nrow = 2, labels = c("(c)", "(d)"))
ACD <- plot_grid(AB[[1]], CD, align = 'h', axis = 'b', ncol = 2, labels = c("(a)", ""), rel_widths = c(2,1))
Fig1 <- plot_grid(ACD, BE, align = 'hv', axis = 'l', nrow = 2, rel_heights = c(3,2))



ggsave("plots/FIG1.pdf", Fig1, width = 13, height = 13)
system("open plots/FIG1.pdf")




