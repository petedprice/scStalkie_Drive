library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
load("data/RData/param_checks/integrated_seurat_nf200_mtr0.20_gu0.RData")
seurat_integrated$treatment <- "SR"
seurat_integrated$treatment[grep("st", seurat_integrated$sample)] <- "ST"
DefaultAssay(seurat_integrated) <- "integrated"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")


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
seurat_integrated <- RunPCA(seurat_integrated, dims = 1:pcs)
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
ggsave("plots/clustering_plots.pdf", clustering_plots, height = 30, width = 30)
ggsave("plots/clustering_plots.png", clustering_plots, height = 30, width = 30)

res <- 'integrated_snn_res.0.4'
Idents(seurat_integrated) <- res
save(seurat_integrated, file = 'data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned.RData')


######### VISUALISING MARKER GENES ------
DefaultAssay(seurat_integrated) <- "RNA"

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
    gsub("_", "-", .) %>% intersect(rownames(seurat_integrated@assays$RNA@counts))
  
  new_name_func <- function(gene){
    nn <- paste(markers$`Drosophila.Marker.(Gene.Name)`[grep(gene, gsub("_", "-", markers$T.dal.Ortholog ))], collapse = ", ")
    nn2 <- paste0(nn , " (", gene, ")")
  }
  
  new_names <- sapply(mk_genes, new_name_func)
  
  dps_list <- lapply(cluster_genes, function(x) DotPlot(seurat_integrated, features = x, assay = "RNA")+coord_flip() +   scale_x_discrete(labels = new_names)
  )
  dps_all <- DotPlot(seurat_integrated, features = unique(unlist(one_cluster_genes)), assay = "RNA")+coord_flip() + 
    scale_x_discrete(labels = new_names)
  dps <- ggarrange(plotlist = dps_list, ncol = 2, nrow = ceiling(length(dps_list)/2), labels = names(cluster_genes))
  #arrange them again but with names as titles
  
  ggsave(paste0("plots/", m, "_tdal_markers_dotplot_clus.pdf"), dps, height = 30, width = 30)
  ggsave(paste0("plots/", m, "_tdal_markers_dotplot_clus.png"), dps, height = 30, width = 30)
  
  ggsave("plots/tdal_markers_dotplot_all_clus.pdf", dps_all, height = 20, width = 14)
  
  dmp <- DimPlot(seurat_integrated, label = T)
  ggsave("plots/tdal_markers_dimplot_clus.pdf", dmp, height = 10, width = 10)
}

#----------------------------------------------


main_figure_markers <- markers$T.dal.Ortholog[markers$Key.Marker == "Yes"] %>% 
  gsub("_", "-", .) %>% 
  gsub(" ", "", .) %>% 
  unique() %>% 
  strsplit(",") %>% unlist() %>% 
  lapply(., function(x)(grep(x, new_names))) %>% unlist() %>% 
  new_names[.] %>% 
  .[c(14,3,1,5,13,4,8,9,2,6,12)]
main_figure_markers <- factor(main_figure_markers, levels = rev(main_figure_markers))
dps_all_figure <- DotPlot(seurat_integrated, features = names(main_figure_markers), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((main_figure_markers)))
DP <- DimPlot(seurat_integrated, label = T)
FP <- FeaturePlot(seurat_integrated, features = 'nFeature_RNA')
FP2 <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.1, y = nFeature_RNA)) + geom_boxplot()
ggarrange(plotlist = list(dps_all_figure, DP, FP, FP2), ncol = 4, nrow = 1, labels = c("A", "B", "C"))


#################### FINAL ASSIGNMENT -------------------------

Unknown <- c(2,4,6,12)
seurat_integrated@meta.data$celltype <- "Keep"

seurat_integrated@meta.data$celltype[which(seurat_integrated@meta.data$integrated_snn_res.0.4 %in% Unknown)] <- "Unknown"
seurat_integrated_ss <- subset(seurat_integrated, idents = c("Keep"))

save(seurat_integrated, file = "data/RData//integrated_seurat_nf200_mtr0.20_gu0_cleaned_celltype.RData")
###############################################################
 


