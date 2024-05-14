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

EBI_markers_out <- read.csv("outdata/markers_EBI.csv") #swap out _ for - in gene
EBI_markers_out$marker <- gsub("_", "-", EBI_markers_out$marker)
  
markers_elife_atlas <- read.csv("outdata/markers_elife_atlas.csv")
markers_elife_atlas$marker <- gsub("_", "-", markers_elife_atlas$marker)

markers_atlaslistdf_out <- read.csv("outdata/markers_atlas.csv")
markers_atlaslistdf_out$marker <- gsub("_", "-", markers_atlaslistdf_out$marker)
markers_elife_out <- read.csv("outdata/markers_elife.csv")
markers_elife_out$marker <- gsub("_", "-", markers_elife_out$marker)
# 
# 
# extra_cells <- filter(EBI_markers_out, celltype %in% c(
#   'accessory_gland', 'adult_fat_cell', 'muscle_cell', 'seminal_vesicle', 'hemocyte'
# )) 
# head(extra_cells)
# #filter extra cells for each celltype, picking top 10 pvals_adj for each 
# extra_cells <- extra_cells %>% 
#   group_by(celltype) %>% 
#   arrange(pvals_adj) %>% 
#   slice(1:15) %>% 
#   ungroup() %>% 
#   distinct()
# 
# elife_plus_extra <- extra_cells %>% 
#   select(celltype, marker) %>%
#   bind_rows(markers_elife_out)
# 
# write.csv(elife_plus_extra, "outdata/markers_elife_plus_extras.csv", quote = F, row.names = F)
# table(elife_plus_extra$celltype)

##################### ---------------------------------------------
  
  
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
    #other_plot <- FeaturePlot(seurat_object, features = genes)
    #ggsave(paste0("plots/markers/", mk_src, "_", c, "_featureplot.pdf"), other_plot, width = 20, height = 20)
  }
  plot <- FeaturePlot(seurat_object, features = plot_features)
  ggsave(paste0("plots/markers/", mk_src, "_featureplot.pdf"), plot, width = 20, height = 20)
}

feature_maker_plotter(EBI_markers_out, seurat_marker, "EBI")
feature_maker_plotter(markers_elife_atlas, seurat_marker, "elife_atlas")
feature_maker_plotter(markers_atlaslistdf_out, seurat_marker, "atlas")
feature_maker_plotter(markers_elife_out, seurat_marker, "elife")
