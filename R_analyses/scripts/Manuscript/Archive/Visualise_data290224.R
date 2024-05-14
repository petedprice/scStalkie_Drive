library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
load("indata/RData/integrated_seurat_nf200_mtr0.20_gu0.RData")
DefaultAssay(seurat_integrated) <- "RNA"


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
pcs
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:pcs)

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
  rename(resolution = integrated_snn_res.) %>%
  ggplot(aes(x = resolution, y = no_clusters)) +
  geom_point()


clustering_plots <- ggarrange(EP, CT, CT_noclusters, ncol = 1, nrow = 3)
ggsave("plots/clustering_plots.pdf", clustering_plots, height = 30, width = 30)

Idents(seurat_integrated) <- "integrated_snn_res.0.4"


######### VISUALISING MARKER GENES ------
markers <- read.xlsx("indata/markers/TDAL_CONSENSUS_MARKERS130224.xlsx", sheet = 2)[,c(1,2,3,4,5,7)] %>% 
  filter(!is.na(Confidence))

clusters <- str_split(markers$Celltype, ",") %>% unlist() %>% #remove space at start of string
  gsub("^ ", "", .) %>% unique()

clusters <- c("Spermatocyte", "Spermatid", "Spermatagonia", "Cyst", "Hub", "GSC")
cluster_genes <- list()
for (c in clusters){
  genes<- markers$T.dal_marker[grepl(tolower(c), tolower(markers$Celltype))]
  genes2 <- genes %>% gsub(",", " ", .) %>% 
    str_split(" ") %>% unlist() %>% gsub("\\(|\\)", "", .) %>% 
    gsub("_", "-", .) %>% intersect(rownames(seurat_integrated@assays$RNA@counts))
  cluster_genes[[c]] <- genes2
  
}

mk_genes <- markers$T.dal_marker %>% gsub(",", " ", .) %>% 
  str_split(" ") %>% unlist() %>% gsub("\\(|\\)", "", .) %>% 
  gsub("_", "-", .) %>% intersect(rownames(seurat_integrated@assays$RNA@counts))


dps_list <- lapply(cluster_genes, function(x) DotPlot(seurat_integrated, features = x, assay = "RNA")+coord_flip())
dps_all <- DotPlot(seurat_integrated, features = unique(unlist(cluster_genes)), assay = "RNA")+coord_flip()
dps <- ggarrange(plotlist = dps1, ncol = 2, nrow = 4, labels = names(cluster_genes))
#arrange them again but with names as titles
ggsave(paste0("plots/tdal_markers_dotplot_clus.pdf"), dp1, height = 20, width = 20)
ggsave(paste0("plots/tdal_markers_dotplot_all_clus.pdf"), dps1_all, height = 20, width = 14)

dmp1 <- DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.1', label = T)
ggsave(paste("plots/tdal_markers_dimplot_clus", dp1_cls_no, ".pdf"), dmp1, height = 10, width = 10)

dp2_cls_no <- names(table(seurat_integrated$integrated_snn_res.0.2)) %>% length()

dps2 <- lapply(cluster_genes, function(x) DotPlot(seurat_integrated, features = x, group.by = 'integrated_snn_res.0.2', assay = "RNA") +coord_flip())
dps2_all <- DotPlot(seurat_integrated, features = unique(unlist(cluster_genes)), group.by = 'integrated_snn_res.0.2', assay = "RNA")+coord_flip()
dp2 <- ggarrange(plotlist = dps2, ncol = 2, nrow = 4, labels = names(cluster_genes))
ggsave(paste0("plots/tdal_markers_dotplot_clus", dp2_cls_no, ".pdf"), dp2, height = 20, width = 20)
ggsave(paste0("plots/tdal_markers_dotplot_all_clus", dp2_cls_no, ".pdf"), dps2_all, height = 20, width = 14)
dmp2 <- DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.2', label = T)
ggsave(paste("plots/tdal_markers_dimplot_clus", dp2_cls_no, ".pdf"), dmp2, height = 10, width = 10)

dp3_cls_no <- names(table(seurat_integrated$integrated_snn_res.0.3)) %>% length()

dps3 <- lapply(cluster_genes, function(x) DotPlot(seurat_integrated, features = x, group.by = 'integrated_snn_res.0.3', assay = "RNA")+coord_flip())
dps3_all <- DotPlot(seurat_integrated, features = unique(unlist(cluster_genes)), group.by = 'integrated_snn_res.0.3', assay = "RNA")+coord_flip()
dp3 <- ggarrange(plotlist = dps3, ncol = 2, nrow = 4, labels = names(cluster_genes))
ggsave(paste0("plots/tdal_markers_dotplot_clus", dp3_cls_no, ".pdf"), dp3, height = 20, width = 20)
ggsave(paste0("plots/tdal_markers_dotplot_all_clus", dp3_cls_no, ".pdf"), dps3_all, height = 20, width = 14)

dmp3 <- DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.3', label = T)
ggsave(paste("plots/tdal_markers_dimplot_clus", dp3_cls_no, ".pdf"), dmp3, height = 10, width = 10)

data_info <- list()
data_info[[1]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.1, fill = treatment, y = mitoRatio)) + 
  geom_boxplot()
data_info[[2]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.1, fill = treatment, y = nFeature_SCT)) + 
  geom_boxplot()
data_info[[3]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.1, fill = treatment, y = nCount_SCT)) + 
  geom_boxplot()
data_info[[4]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.1, fill = treatment, y = xprop)) + 
  geom_boxplot()
ggpubr::ggarrange(plotlist = data_info, ncol = 2, nrow = 2, labels = c(1:4))
ggsave("plots/expression_info.pdf", height = 20, width = 30)


seurat_integrated@meta.data %>% 
  ggplot(aes(x = xprop, y = nFeature_SCT)) +
  geom_smooth(method = 'lm') + ggpubr::stat_cor() +
  geom_point()
DefaultAssay(seurat_integrated) <- "SCT"
seurat_integrated <- AverageExpression(seurat_integrated)

HtM <- DoHeatmap(seurat_integrated, features = unique(unlist(cluster_genes)), group.by = 'celltype')
ggsave("plots/tdal_markers_heatmap.pdf", height = 20, width = 20)

ortholog_table$status = "No Ortholog" #add new variable where if OMA_GENE is not NA you give it "OMA", 
#if OF_DMEL is not NA you give it "OF", if both you give it "both"
ortholog_table$status[which(!is.na(ortholog_table$OMA_REFGENE) & is.na(ortholog_table$OF_DMEL))] = "OMA"
ortholog_table$status[which(!is.na(ortholog_table$OF_DMEL) & is.na(ortholog_table$OMA_REFGENE))] = "OF"
ortholog_table$status[which(!is.na(ortholog_table$OMA_REFGENE) & !is.na(ortholog_table$OF_DMEL))] = "both"


ortho_plot <- ortholog_table %>% select(OMA_REFGENE, OF_DMEL, REF_GENE_NAME) %>% unique()
ortho_plot$status = "No Ortholog" #add new variable where if OMA_GENE is not NA you give it "OMA",
#if OF_DMEL is not NA you give it "OF", if both you give it "both"
ortho_plot$status[which(!is.na(ortho_plot$OMA_REFGENE) & is.na(ortho_plot$OF_DMEL))] = "OMA"
ortho_plot$status[which(!is.na(ortho_plot$OF_DMEL) & is.na(ortho_plot$OMA_REFGENE))] = "Orthofinder"
ortho_plot$status[which(!is.na(ortho_plot$OMA_REFGENE) & !is.na(ortho_plot$OF_DMEL))] = "Both"


ortho_plot %>% ggplot(aes(x = status, fill = status)) + geom_bar() + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Orthologs in TDAL markers", x = "Ortholog Status", y = "Count") +#add colour blind friendly colour scheme
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
ggsave("plots/ortholog_status.pdf", height = 7, width = 7)


#using the cluster_genes list make a data.frame with a column called cluster and a second called gene
# the cluster is the name of the element in the list and the gene will be the different genes under that name
cluster_genes_df <- data.frame(cluster = rep(names(cluster_genes), sapply(cluster_genes, length)), 
                               gene = unlist(cluster_genes))
cluster_genes_df %>% 
  ggplot(aes(x = cluster, fill = cluster)) + geom_bar() +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Number of genes in each cluster", x = "Cluster", y = "Count") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
ggsave("plots/cluster_gene_count.pdf", height = 7, width = 7)

####### del 
DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.1', label = T)
ggplot(seurat_integrated@meta.data, aes(x = integrated_snn_res.0.1, fill = treatment)) + geom_bar(position = 'dodge')
