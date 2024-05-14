library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(ggpubr)
#library(scater)
library(scuttle)
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(data.table)
library("biomaRt")
library(clusterProfiler)
library(ggrepel)
library(lme4)
library(GenomicFeatures)

############### LOAD DATA ----------------------
load("data/RData/marker_seurat.RData")
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
seurat_marker$celltype = seurat_marker$sctype_labels
seurat_marker$treatment = "ST" #which samples start with SR treatment is SR
seurat_marker$treatment[grepl("sr", seurat_marker$sample)] = "SR"
ortholog_table$length <- ortholog_table$end - ortholog_table$start
ortholog_table$mid <- ortholog_table$start + (ortholog_table$length/2)


pg_markers <- readxl::read_excel("data/PopGroup_markers.xlsx") 
colnames(pg_markers)[1] <- "seurat_markers"

for (i in 1:nrow(pg_markers)){
  seurat_marker$celltype[seurat_marker$seurat_clusters %in% pg_markers$seurat_markers[i]] <- pg_markers$Consensus[i]
}
seurat_marker <- subset(seurat_marker, celltype != 'Unknown' & celltype != "Accessory Gland?" & celltype != "Muscle")
seurat_marker$celltype[seurat_marker$`integrated_snn_res.1.5` == 32] <- "Spermatids"
PG_UMAP <- DimPlot(seurat_marker,  label = TRUE, group.by = 'celltype', split.by = 'treatment')
ggsave("plots/UMAP_PopGroup.pdf", PG_UMAP, width = 16, height = 10)
seurat_marker$celltype <- gsub(" ", "_", seurat_marker$celltype)
seurat_marker$celltype <- gsub("/", "_and_", seurat_marker$celltype)
seurat_marker$celltype <- gsub("muscle", "Muscle_cells", seurat_marker$celltype)



seurat_marker$comp_cluster <- seurat_marker$integrated_snn_res.0.3
DimPlot(seurat_marker, group.by = 'comp_cluster', label = T)


EBI_markers_out <- read.csv("outdata/markers_EBI.csv") #swap out _ for - in gene
EBI_markers_out$marker <- gsub("_", "-", EBI_markers_out$marker)

markers_elife_atlas <- read.csv("outdata/markers_elife_atlas.csv")
markers_elife_atlas$marker <- gsub("_", "-", markers_elife_atlas$marker)

markers_atlaslistdf_out <- read.csv("outdata/markers_atlas.csv")
markers_atlaslistdf_out$marker <- gsub("_", "-", markers_atlaslistdf_out$marker)
markers_elife_out <- read.csv("outdata/markers_elife.csv")
markers_elife_out$marker <- gsub("_", "-", markers_elife_out$marker)

wei_markers <- read.table("indata/markers/wei2022_markerlist.txt") %>% 
  merge(ortholog_table, by.x = "V1", by.y = "consensus_gene", all.x = T) %>% 
  mutate(REF_GENE_NAME = gsub("_", "-", REF_GENE_NAME))

seurat_marker_rename <- seurat_marker
for (i in 1:nrow(seurat_marker@assays$RNA)){
  rnm <- rownames(seurat_marker_rename@assays$RNA@data)[i] %>% gsub("-", "_", .)
  chan <- ortholog_table$consensus_gene[which(ortholog_table$REF_GENE_NAME == rnm)][1] %>% tolower() 
  if (is.na(chan) == F){
    rownames(seurat_marker_rename@assays$RNA@data)[i] <- chan
  }
}


DefaultAssay(seurat_marker_rename) <- "RNA"
gs <- intersect(rownames(seurat_marker_rename), tolower(unique(wei_markers$V1))) 
filter(ortholog_table, tolower(consensus_gene) %in% gs) %>% 
  dplyr::select(consensus_gene) %>% unique() %>% dim()


dmp <- DimPlot(seurat_marker_rename, group.by = 'celltype', label = T)
dp <- DotPlot(seurat_marker_rename, group.by = 'celltype', features = tolower(unique(wei_markers$V1)))
wei_celltypes <- ggarrange(dmp, dp, ncol = 1, nrow = 2)
ggsave("plots/Cell_types/wei_celltypes.pdf", wei_celltypes, width = 10, height = 20)
wei_fp <- FeaturePlot(seurat_marker_rename, features = tolower(unique(wei_markers$V1)), ncol = 4)
ggsave("plots/Cell_types/wei_fp.pdf", wei_fp, width = 20, height = 20)
dmp2 <- DimPlot(seurat_marker_rename, group.by = 'comp_cluster', label = T)
dp2 <- DotPlot(seurat_marker_rename, group.by = 'comp_cluster', features = tolower(unique(wei_markers$V1)))
wei_clusters <- ggarrange(dmp2, dp2, ncol = 1, nrow = 2)
ggsave("plots/Cell_types/wei_clusters.pdf", wei_clusters, width = 10, height = 20)


########### WITT MARKERS -------------
witt_markers <- read.table("indata/markers/witt2019_markerlist.txt") %>% 
  merge(ortholog_table, by.x = "V1", by.y = "consensus_gene", all.x = T) %>% 
  mutate(REF_GENE_NAME = gsub("_", "-", REF_GENE_NAME))
gs <- intersect(rownames(seurat_marker_rename), tolower(unique(witt_markers$V1)))
dmp <- DimPlot(seurat_marker_rename, group.by = 'celltype', label = T)
dp <- DotPlot(seurat_marker_rename, group.by = 'celltype', features = tolower(unique(witt_markers$V1)))
witt_celltypes <- ggarrange(dmp, dp, ncol = 1, nrow = 2)
ggsave("plots/Cell_types/witt_celltypes.pdf", witt_celltypes, width = 10, height = 20)
witt_fp <- FeaturePlot(seurat_marker_rename, features = tolower(unique(witt_markers$V1)), ncol = 2)
ggsave("plots/Cell_types/witt_fp.pdf", witt_fp, width = 20, height = 20)
dmp2 <- DimPlot(seurat_marker_rename, group.by = 'comp_cluster', label = T)
dp2 <- DotPlot(seurat_marker_rename, group.by = 'comp_cluster', features = tolower(unique(witt_markers$V1)))
witt_clusters <- ggarrange(dmp2, dp2, ncol = 1, nrow = 2)
ggsave("plots/Cell_types/witt_clusters.pdf", witt_clusters, width = 10, height = 20)


#### MAH MARKERS 
mah_markers <- read.table("outdata/markers_mah.csv", sep = ',', header = T) %>% 
  merge(ortholog_table, by.x = "marker", by.y = "REF_GENE_NAME", all.x = T) 

dmp <- DimPlot(seurat_marker_rename, group.by = 'celltype', label = T)
dp <- DotPlot(seurat_marker_rename, group.by = 'celltype', features = unique(mah_markers$consensus_gene))
mah_celltypes <- ggarrange(dmp, dp, ncol = 1, nrow = 2)
ggsave("plots/Cell_types/mah_celltypes.pdf", mah_celltypes, width = 20, height = 30)

mah_fp <- FeaturePlot(seurat_marker_rename, features = mah_markers$consensus_gene, ncol = 2)
ggsave("plots/Cell_types/mah_fp.pdf", mah_fp, width = 20, height = 20)

dmp2 <- DimPlot(seurat_marker_rename, group.by = 'comp_cluster', label = T)
dp2 <- DotPlot(seurat_marker_rename, group.by = 'comp_cluster', features = unique(mah_markers$consensus_gene))
mah_clusters <- ggarrange(dmp2, dp2, ncol = 1, nrow = 2)
ggsave("plots/Cell_types/mah_clusters.pdf", mah_clusters, width = 20, height = 30)

#########








############ GET MARKERS/CLUSTERS ---------------------

UMAP <- DimPlot(seurat_marker, group.by = 'seurat_clusters', label = T)
ggsave('plots/UMAP.pdf', UMAP)

features = markers_atlaslistdf_out
seurat_object = seurat_marker
mk_src = "atlaslistdf"
feature_maker_plotter <- function(features, seurat_object, mk_src){
  DefaultAssay(seurat_object) <- "RNA"
  plot_features <- c()
  for (c in unique(features$celltype)){
    filtf <- filter(features, celltype == c)
    genes = intersect(filtf$marker, rownames(seurat_object))
    if (length(genes) == 0){
      next
    }
    mexp <- PercentageFeatureSet(object = seurat_object, features = genes)
    fname = paste("mexp_", c, sep = "")
    seurat_object <- AddMetaData(object = seurat_object, metadata = mexp, col.name = fname)
    plot_features <- c(plot_features, fname)
    other_plot <- FeaturePlot(seurat_object, features = genes)
    ggsave(paste0("plots/markers/", mk_src, "_", c, "_featureplot.pdf"), other_plot, width = 20, height = 20)
    op2 <- DotPlot(seurat_object, features = genes, group.by = 'comp_cluster')
    ggsave(paste0("plots/markers/", mk_src, "_", c, "_dotplot.pdf"), op2, width = 20, height = 20)
  }
  plot <- FeaturePlot(seurat_object, features = plot_features)
  ggsave(paste0("plots/markers/", mk_src, "_featureplot.pdf"), plot, width = 20, height = 20)
}

feature_maker_plotter(EBI_markers_out, seurat_marker, "EBI")
feature_maker_plotter(markers_elife_atlas, seurat_marker, "elife_atlas")
feature_maker_plotter(markers_atlaslistdf_out, seurat_marker, "atlas")
feature_maker_plotter(markers_elife_out, seurat_marker, "elife")
feature_maker_plotter(mah_markers, seurat_marker, "mah")

