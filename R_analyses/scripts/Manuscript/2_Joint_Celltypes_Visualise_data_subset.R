library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
rm(list = ls())


ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

recluster=FALSE

if (recluster == TRUE){
  load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_subset.RData")
  
  ######################## RE Run but REMOVE CLUSTERS -----------
  
  DefaultAssay(seurat_integrated_ss) <- "integrated"
  
  seurat_integrated_ss <- RunPCA(seurat_integrated_ss, dims = 1:40)
  seurat_integrated_ss <- RunUMAP(seurat_integrated_ss, dims = 1:40)
  
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
  
  
  save(seurat_integrated_ss, file = 'data/RData/integrated_seurat_nf200_mtr0.20_gu0_reclus_subset.RData')
 } else {
  load('data/RData/integrated_seurat_nf200_mtr0.20_gu0_reclus_subset.RData')
 }

res <- 'integrated_snn_res.0.4'
Idents(seurat_integrated_ss) <- res
######### VISUALISING MARKER GENES ------
DefaultAssay(seurat_integrated_ss) <- "RNA"

markers <- read.xlsx("data/MANUSCRIPT/Drive Manuscript Supplementary Tables.xlsx", sheet = 2)[,c(1:6)] %>% 
  filter(`DOI.(Source)` != "https://doi.org/10.1038/s41467-021-20897-y" & !is.na(Cell.type))

m = "Cell.type"
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
    gsub("_", "-", .) %>% intersect(rownames(seurat_integrated_ss@assays$RNA@counts))
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


Figure_marker_names_df <- data.frame(matrix(nrow = 33))
Figure_marker_names_df$tdal <- markers$T.dal.Ortholog %>% 
  gsub("_", "-", .) %>% 
  gsub(" ", "", .) %>% 
  unique() %>% 
  strsplit(",") %>% unlist()
Figure_marker_names_df$dros <- Figure_marker_names_df$tdal %>% 
  lapply(., function(x)(grep(x, new_names))) %>% unlist() %>% 
  new_names[.]

main_figure_markers_indx <- c(20,25,24,23,22,21,10,13,29,17,27,28,31)
main_figure_markers <- Figure_marker_names_df$dros[main_figure_markers_indx]
names(main_figure_markers) <- Figure_marker_names_df$tdal[main_figure_markers_indx]
main_figure_markers <- factor(main_figure_markers, levels = rev(main_figure_markers))


all_figure_markers_indx <- c(20,33,25,24,23,22,21,2,3,4,7,8,10,11,12,13,15,30,17,26,27,14,28,31)
all_figure_markers <- Figure_marker_names_df$dros[all_figure_markers_indx]
names(all_figure_markers) <- Figure_marker_names_df$tdal[all_figure_markers_indx]

### Plotting key markers against numbered cell clusters for assignment ---
Idents(seurat_integrated_ss) <- seurat_integrated_ss$integrated_snn_res.0.4
dps_all_figure <- DotPlot(seurat_integrated_ss, features = names(main_figure_markers), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((main_figure_markers)))
nfeatures_figure <- seurat_integrated_ss@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = nFeature_RNA)) +
  geom_boxplot() + 
  theme_bw() + scale_fill_discrete(name = "Celltype")
UMAP <- DimPlot(seurat_integrated_ss, label = T)
UMAP2 <- DimPlot(seurat_integrated_ss, group.by = "Phase")
ggarrange(dps_all_figure, nfeatures_figure, UMAP, UMAP2, ncol = 4) %>% 
  ggsave( "plots/subset_celltype_key_markers.pdf", ., width = 40, height = 10)
#--------------------------------------------------------------------------

#0.4 clusters
if (res == "integrated_snn_res.0.4"){
  Muscle <- c(15)
  Pre_meiotic_cyst <- c(6,7,11,0)
  Post_meiotic_cyst <- c(4,10)
  `GSC/Spermatogonia` <- c(1)
  Primary_spermatocytes <- c(3,14)
  Secondary_spermatocytes <- c(12,13)
  Spermatids <- c(2,5,8,9)
} else if (res == "integrated_snn_res.0.75"){
  #0.75 clusters
  Muscle <- c(18)
  Pre_meiotic_cyst <- c(6,7,0,4)
  Post_meiotic_cyst <- c(19,17,2,11)
  `GSC/Spermatogonia` <- c(1,3,12)
  Primary_spermatocytes <- c(20,22) ### FINISH ASSIGNMENT ##
  Secondary_spermatocytes <- c(9,16)
  Spermatids <- c(22,15,3,10,11)
}

#seurat_final <- CellSelector(DimPlot(seurat_final), seurat_final, "outliers1")
#seurat_final <- CellSelector(DimPlot(seurat_final), seurat_final, "outliers2")


#DotPlot(seurat_final, features = names(main_figure_markers), assay = "RNA")+coord_flip() +
#  scale_x_discrete(labels = as.vector((main_figure_markers)))

seurat_integrated_ss@meta.data$celltype <- 'NA'
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Muscle] <- "Muscle"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Primary_spermatocytes] <- "Primary spermatocytes"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Secondary_spermatocytes] <- "Secondary spermatocytes"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatids] <- "Spermatids"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% `GSC/Spermatogonia`] <- "GSC/Spermatogonia"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Pre_meiotic_cyst] <- "Early cyst"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Post_meiotic_cyst] <- "Late cyst"


seurat_final <- seurat_integrated_ss

Idents(seurat_final) <- "celltype"
levels(seurat_final) <- c("Muscle", "Early cyst", "Late cyst", 
                                  "GSC/Spermatogonia", "Primary spermatocytes", 
                                  "Secondary spermatocytes", "Spermatids")


save(seurat_final, main_figure_markers, all_figure_markers, file = 'data/RData/seurat_final.RData')



