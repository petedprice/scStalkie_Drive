library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
library(cowplot)
rm(list = ls())

load("data/RData/integrated_seurat_nf200_mtr1_gu0.RData")
seurat_integrated$treatment <- "SR"
seurat_integrated$treatment[grep("st", seurat_integrated$sample)] <- "ST"
DefaultAssay(seurat_integrated) <- "integrated"
ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

seurat_integrated@meta.data %>% 
  mutate(mitoclass = ifelse(mitoRatio > 0.2, "high mitochondrial expression", "pass")) %>% 
  group_by(treatment, mitoclass) %>% 
  summarise(n = dplyr::n())

seurat_integrated@meta.data %>% 
  mutate(mitoclass = ifelse(mitoRatio > 0.2, 1, 0)) %>% 
  lme4::glmer(mitoclass ~ treatment + (1|sample), data =., family = 'binomial') %>% 
  summary()



#############################################

seurat_integrated <- RunPCA(seurat_integrated)

######### Determining number of PCAs to use for clustering/UMAP etc.
EP <- ElbowPlot(seurat_integrated, ndims = 40)
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
pcs <- min(co1, co2)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:pcs)
seurat_integrated <- RunTSNE(seurat_integrated, dims = 1:pcs)

######### Determining resolution to use --------
seurat_integrated <- FindNeighbors(object=seurat_integrated, dims=1:pcs)
seurat_integrated <- FindClusters(object=seurat_integrated, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                                                           0.75, 1, 1.25, 1.5, 1.75, 2, 
                                                                           2.5, 3))
clusttree <- clustree(seurat_integrated)
CT <- plot(clusttree)
CT_noclusters <- clusttree$data %>% 
  group_by(integrated_snn_res.) %>% 
  summarise(no_clusters = n_distinct(cluster)) %>% 
  dplyr::rename(resolution = `integrated_snn_res.`) %>%
  ggplot(aes(x = resolution, y = no_clusters)) +
  geom_point()


clustering_plots <- ggarrange(EP, CT, CT_noclusters, ncol = 1, nrow = 3)

res <- 'integrated_snn_res.0.4'
Idents(seurat_integrated) <- res




################ MITO EXPRESSION 



seurat_integrated@meta.data$mitoRatio <- seurat_integrated@meta.data$mitoRatio *100
a <- FeaturePlot(subset(seurat_integrated, treatment == "ST"), features = 'mitoRatio')
b <- FeaturePlot(subset(seurat_integrated, treatment == "SR"), features = 'mitoRatio')

c <- FeaturePlot(subset(seurat_integrated, treatment == "ST"), features = 'nFeature_RNA')
d <- FeaturePlot(subset(seurat_integrated, treatment == "SR"), features = 'nFeature_RNA')


label1 <- labs(x = "UMAP 1", y = "UMAP 2", color = "Mitochondrial \nexpression (%)")
label2 <- labs(x = "UMAP 1", y = "UMAP 2", color = "Number of \nfeatures")

a2  <- a+ label1 + labs(title = "ST") 
b2 <- b + label1 + labs(title = "SR") 
c2 <- c + label2 + labs(title = "ST")
d2 <- d + label2 + labs(title = "SR")

mitoplot1 <- ggarrange(a2,b2, labels = c("(a)", "(b)"), common.legend = T, legend = 'right')
#ggsave("plots/manuscript_plots/S12.tiff", height = 5, width = 9.5)
#system("open plots/manuscript_plots/S12.tiff")

mitoplot2 <- ggarrange(c2,d2, labels = c("(c)", "(d)"), common.legend = T, legend = 'right')
mitoplotall <- ggarrange(mitoplot1, mitoplot2, nrow = 2)
ggsave("plots/manuscript_plots/S12.tiff", height = 10, width = 9.5)
system("open plots/manuscript_plots/S12.tiff")



####################################################################################################################

########################## EXTRA CHECKS ########################## 

####################################################################################################################

######### VISUALISING MARKER GENES ------
markers <- read.xlsx("data/MANUSCRIPT/Drive Manuscript Supplementary Tables.xlsx", sheet = 2)[1:31,c(1:6)] 
row_order <- c("M", "CYSC", "EC", "LC", "C", "E", "G", "G, PS", "PS, SS", "PS, SS, ST", "PS, ST", "ST", "LST")
markers$Drosophila.cell.type <- factor(markers$Drosophila.cell.type, levels = row_order)
markers <- markers[order(markers$Drosophila.cell.type),]

mk_genes <- markers$Teleopsis.dalmanni.Ortholog %>% gsub(",", " ", .) %>% 
  str_split(" ") %>% unlist() %>% 
  gsub("_", "-", .) %>% intersect(rownames(seurat_integrated@assays$RNA))


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

main_figure_markers <- all_figure_markers[c(1,3:7,14,15,18,21,23,25,24,26,28,29)]


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
dps_all_figure <- DotPlot(seurat_integrated, features = names(mfms), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((mfms))) + 
  labs(y = '', x = "Marker gene expression") + 
  theme_classic() + scale_colour_gradientn(colours = colorspace::diverge_hcl(8)) + 
  theme(legend.position="right", legend.box = "vertical", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) + 
  guides(color = guide_colorbar(title = 'Scaled expression'), 
         size = guide_legend(title = '% of cells expressed'))


nfeatures_figure <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = nFeature_RNA, fill = treatment)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + 
  labs(x = "", y = "N\u00b0 expressed  \n genes") + 
  theme(legend.position="none")


UMAP <- DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.4', label = T) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(legend.position="none")
UMAP_final <- UMAP[[1]]$data %>% 
  ggplot(aes(x = umap_1, y = umap_2, colour = integrated_snn_res.0.4)) + 
  geom_point(stroke = NA, size = .7) +  
  guides(colour = guide_legend(override.aes = list(size=3), title = "Cell type")) + 
  theme_classic() + labs(x = "UMAP 1", y = "UMAP 2", legend = 'Cell type') +
  theme(legend.position = 'none', legend.box = "vertical", 
        axis.text.x = element_text(color="black"), 
        axis.ticks = element_line(color = "black")) +
  shadowtext::geom_shadowtext(data = UMAP[[1]]$layers[[2]]$data, aes(label = integrated_snn_res.0.4), colour = 'black', 
                              size = 4, nudge_y = -0.5, nudge_x = 0, 
                              bg.color = 'white', bg.r = 0.1, bg.size = 0.5, 
                              fontface = "bold") 


mitoExpression_fig  <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = mitoRatio * 100, fill = treatment)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() +
  labs(x = "", y = "% mitochondrial expression") + 
  theme(legend.position="top")

FP <- FeaturePlot(seurat_integrated, features = 'mitoRatio', split.by = 'treatment')

mito_exp_figure <- ggarrange(ggarrange(UMAP_final, mitoExpression_fig, nfeatures_figure, dps_all_figure), 
          FP, nrow = 2, heights = c(2,1))

ggsave("plots/no_mito_filter_summary_plot.pdf", height = 15, width = 14)
system('open plots/no_mito_filter_summary_plot.pdf')

