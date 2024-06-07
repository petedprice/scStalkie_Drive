library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")


load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned_celltype.RData")

######################## RE Run but REMOVE CLUSTERS -----------


Idents(seurat_integrated) <- seurat_integrated$celltype
seurat_integrated_ss <- subset(seurat_integrated, idents = c("TBC", "Cyst", "GSC/Spermatogonia", "Spermatid", "Spermatocyte", "Muscle"))
DefaultAssay(seurat_integrated_ss) <- "integrated"

######### Determining number of PCAs to use for clustering/UMAP etc.
EP <- ElbowPlot(seurat_integrated_ss, ndims = 40)
pct <- seurat_integrated_ss[["pca"]]@stdev / sum(seurat_integrated_ss[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
pcs <- min(co1, co2)
pcs
seurat_integrated_ss <- RunUMAP(seurat_integrated_ss, dims = 1:pcs)
seurat_integrated_ss <- RunPCA(seurat_integrated_ss, dims = 1:pcs)


######### Determining resolution to use --------
seurat_integrated_ss <- FindNeighbors(object=seurat_integrated_ss, dims=1:pcs)
seurat_integrated_ss <- FindClusters(object=seurat_integrated_ss, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                                                                 0.75, 1, 1.25, 1.5, 1.75, 2, 
                                                                                 2.5, 3))
clusttree <- clustree(seurat_integrated_ss)
CT <- plot(clusttree)
CT_noclusters <- clusttree$data %>% 
  group_by(integrated_snn_res.) %>% 
  summarise(no_clusters = n_distinct(cluster)) %>% 
  rename(resolution = integrated_snn_res.) %>%
  ggplot(aes(x = resolution, y = no_clusters)) +
  geom_point()


clustering_plots <- ggarrange(EP, CT, CT_noclusters, ncol = 1, nrow = 3)
ggsave("plots/clustering_plots_ss.pdf", clustering_plots, height = 30, width = 30)
ggsave("plots/clustering_plots_ss.png", clustering_plots, height = 30, width = 30)

res <- 'integrated_snn_res.0.4'
Idents(seurat_integrated_ss) <- res

save(seurat_integrated_ss, file = 'data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned_ss.RData')


######### VISUALISING MARKER GENES ------
DefaultAssay(seurat_integrated_ss) <- "RNA"

markers <- read.xlsx("data/MANUSCRIPT/DRIVE_MANUSCRIPT_SUPPLEMENTARY_TABLES.xlsx", sheet = 1)[,c(2,3,4,5,6,7)] %>% 
  filter(`DOI.(Source)` != "https://doi.org/10.1038/s41467-021-20897-y")

for (m in c("Celltype.(Specific)", "Celltype.(Broad)")){
  clusters <- str_split(markers[,m], ",") %>% unlist() %>% #remove space at start of string
    gsub("^ ", "", .) %>% unique()
  clusters <- clusters[!is.na(clusters)]
  #clusters <- markers$Medium_Group %>% unique()
  clusters <- markers[,m] %>% unique()
  
  clusters <- clusters[!is.na(clusters)]
  cluster_genes <- list()
  
  for (c in clusters){
    genes <- markers$T.dal.Ortholog[markers[,m] == c]
    #genes<- markers$T.dal_marker[grepl(tolower(c), tolower(markers$Consensus_Marker))]
    genes2 <- genes %>% gsub(",", " ", .) %>% 
      str_split(" ") %>% unlist() %>% gsub("\\(|\\)", "", .) %>% 
      gsub("_", "-", .) %>% intersect(rownames(seurat_integrated@assays$RNA@counts))
    cluster_genes[[c]] <- genes2
    
  }
  cluster_genes <- compact(cluster_genes)
  one_cluster_genes_func <- function(cluster, cluster_genes){
    other_clus <- names(cluster_genes) %>% setdiff(cluster)
    all_genes <- unlist(cluster_genes[other_clus])
    one_clus <- cluster_genes[[cluster]][which(!(cluster_genes[[cluster]] %in% all_genes))]
    if (length(one_clus) > 0){
      return(one_clus)
    }
  }
  
  one_cluster_genes <- sapply(clusters, one_cluster_genes_func, cluster_genes = cluster_genes, 
                              simplify = F, USE.NAMES = T) %>% 
    compact()
  
  ##### OVERLAP GENES ----
  
  
  
  mk_genes <- markers$T.dal.Ortholog %>% gsub(",", " ", .) %>% 
    str_split(" ") %>% unlist() %>% gsub("\\(|\\)", "", .) %>% 
    gsub("_", "-", .) %>% intersect(rownames(seurat_integrated_ss@assays$RNA@counts))
  
  new_name_func <- function(gene){
    nn <- paste(markers$`Drosophila.Marker.(Gene.Name)`[grep(gene, gsub("_", "-", markers$T.dal.Ortholog ))], collapse = ", ")
    nn2 <- paste0(nn , " (", gene, ")")
  }
  
  new_names <- sapply(mk_genes, new_name_func)
  
  dps_list <- lapply(cluster_genes, function(x) DotPlot(seurat_integrated_ss, features = x, assay = "RNA")+coord_flip() +   scale_x_discrete(labels = new_names)
  )
  dps_all <- DotPlot(seurat_integrated_ss, features = unique(unlist(one_cluster_genes)), assay = "RNA")+coord_flip() + 
    scale_x_discrete(labels = new_names)
  dps <- ggarrange(plotlist = dps_list, ncol = 2, nrow = ceiling(length(dps_list)/2), labels = names(cluster_genes))
  #arrange them again but with names as titles
  
  ggsave(paste0("plots/", m, "_tdal_markers_dotplot_clus_ss.pdf"), dps, height = 30, width = 30)
  ggsave(paste0("plots/", m, "_tdal_markers_dotplot_clus_ss.png"), dps, height = 30, width = 30)
  
  ggsave("plots/tdal_markers_dotplot_all_clus_ss.pdf", dps_all, height = 20, width = 14)
  
  dmp <- DimPlot(seurat_integrated_ss, label = T)
  ggsave("plots/tdal_markers_dimplot_clus_ss.pdf", dmp, height = 10, width = 10)
}

main_figure_markers <- markers$T.dal.Ortholog[markers$Key.Marker == "Yes"] %>% 
  gsub("_", "-", .) %>% 
  gsub(" ", "", .) %>% 
  unique() %>% 
  strsplit(",") %>% unlist() %>% 
  lapply(., function(x)(grep(x, new_names))) %>% unlist() %>% 
  new_names[.] %>% 
  .[c(10,1,5,9,3,4,6,2,8)]
main_figure_markers <- factor(main_figure_markers, levels = rev(main_figure_markers))

Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype

### Plotting key markers against numbered cell clusters for assignment ---

Idents(seurat_integrated_ss) <- seurat_integrated_ss$integrated_snn_res.0.4
dps_all_figure <- DotPlot(seurat_integrated_ss, features = names(main_figure_markers), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((main_figure_markers)))
nfeatures_figure <- seurat_integrated_ss@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = nFeature_RNA, fill = integrated_snn_res.0.4)) +
  geom_boxplot() + 
  theme_bw() + scale_fill_discrete(name = "Celltype")
UMAP <- DimPlot(seurat_integrated_ss, label = T)
ggarrange(dps_all_figure, UMAP, ncol = 2) %>% 
  ggsave( "plots/subset_celltype_key_markers.pdf", ., width = 20, height = 10)

#--------------------------------------------------------------------------

Muscle <- c(17)
Spermatocytes <- c(10,11,13)
Spermatids <- c(3,4)
`GSC/Spermatogonia` <- c(0,2,7)
Pre_meiotic_cyst <- c(5,12,14,15)
Post_meiotic_cyst <- c(1,6,8,9,16)

### Plotting key markers against numbered cell clusters for assignment ---

seurat_integrated_ss@meta.data$celltype <- 'NA'
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Muscle] <- "Muscle"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatocytes] <- "Spermatocytes"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatids] <- "Spermatids"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% `GSC/Spermatogonia`] <- "GSC/Spermatogonia"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Pre_meiotic_cyst] <- "Pre-meiotic \ncyst"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Post_meiotic_cyst] <- "Post-meiotic \ncyst"

Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype
dps_all_figure <- DotPlot(seurat_integrated_ss, features = names(main_figure_markers), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((main_figure_markers)))
dps_all_figure
nfeatures_figure <- seurat_integrated_ss@meta.data %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA, fill = celltype)) +
  geom_boxplot() + 
  theme_bw() + scale_fill_discrete(name = "Celltype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
UMAP <- DimPlot(seurat_integrated_ss, label = T) + 
  labs(x = "UMAP 1", y = "UMAP 2")


#Arrange UMAP, Dotplot and nfeatures_figure with UMAP as first row and twice the height of the other two plots 
ggarrange(UMAP, ggarrange(dps_all_figure, nfeatures_figure, ncol = 2, labels = c("B", "C")), ncol = 1, nrow = 2, heights = c(3,3), labels = c("A")) %>% 
  ggsave( "plots/subset_celltype_key_markers.pdf", ., width = 25, height = 20)

ggarrange(UMAP, dps_all_figure, ncol = 2, nrow = 1, widths = c(4,5)) %>% 
  ggsave( "plots/subset_celltype_key_markers.pdf", ., width = 25, height = 10)



######### ADDITIONAL CELLTYPE ASSIGNMENT USING NO GENES EXPRESSED -----
no_genes_plots <- list()
no_genes_plots[[1]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = nFeature_RNA, fill = treatment)) +
  geom_boxplot()
no_genes_plots[[2]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = log(nCount_RNA), fill = treatment)) +
  geom_boxplot()
no_genes_plots[[3]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = log10GenesPerUMI, fill = treatment)) +
  geom_boxplot()
no_genes_plots[[4]] <- DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.4', label = T)
no_genes_plots_arranged <- ggarrange(plotlist = no_genes_plots, ncol = 2, nrow = 2)
no_genes_plots_arranged
ggsave("plots/no_genes_plots.png", no_genes_plots_arranged, height = 13, width = 15)

no_genes_featplots <- FeaturePlot(seurat_integrated, features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI"), pt.size = 0.5, split.by = 'treatment')
no_genes_featplots <- no_genes_featplots & theme(legend.position = c(0.1,0.2))
ggsave("plots/no_genes_featplots.png", no_genes_featplots, height = 10, width = 15)

######### ADDITIONAL CELLTYPE ASSIGNMENT USING X GENES EXPRESSED -----
Xgenes <- filter(ortholog_table, chr == "Chr_X") %>% 
  dplyr::select(REF_GENE_NAME) %>% unlist() %>% 
  gsub("_", "-", .) %>% 
  c() %>% intersect(rownames(seurat_integrated@assays$RNA@counts))

Xexp <- PercentageFeatureSet(seurat_integrated, features = Xgenes)
seurat_integrated <- AddMetaData(object = seurat_integrated, metadata = Xexp, col.name = 'Xexp')

Xgene_prop <- colSums(seurat_integrated@assays$RNA$counts[Xgenes,] > 0)/
  colSums(seurat_integrated@assays$RNA$counts > 0)

seurat_integrated <- AddMetaData(object = seurat_integrated, metadata = Xgene_prop, col.name = 'Xprop')

xplots <- list()
xplots[[1]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = Xexp, fill = treatment)) +
  geom_boxplot()
xplots[[2]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = Xprop, fill = treatment)) +
  geom_boxplot()
xplots[[3]] <- DimPlot(seurat_integrated, group.by = res, label = T)
xplots[[4]] <- FeaturePlot(seurat_integrated, features = c("Xexp", "Xprop"), pt.size = 0.5, split.by = 'treatment')
xplots_arranged <- ggarrange(plotlist = xplots, ncol = 2, nrow = 2)
xplots_arranged
ggsave("plots/xplots.png", xplots_arranged, height = 13, width = 15)
###

data_info <- list()
data_info[[1]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, fill = treatment, y = mitoRatio)) + 
  geom_boxplot()
data_info[[2]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, fill = treatment, y = nFeature_SCT)) + 
  geom_boxplot()
data_info[[3]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, fill = treatment, y = nCount_SCT)) + 
  geom_boxplot()
data_info[[4]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, fill = treatment, y = xprop)) + 
  geom_boxplot()
ggpubr::ggarrange(plotlist = data_info, ncol = 2, nrow = 2, labels = c(1:4))
ggsave("plots/expression_info.pdf", height = 20, width = 30)





