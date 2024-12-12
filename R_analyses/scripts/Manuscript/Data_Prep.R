library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
rm(list = ls())

load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_newref_Seurat5.1_df0.08.RData")

seurat_integrated$treatment <- "SR"
seurat_integrated$treatment[grep("st", seurat_integrated$sample)] <- "ST"
DefaultAssay(seurat_integrated) <- "integrated"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)



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
ggsave("plots/clustering_plots.png", clustering_plots, height = 30, width = 30)

res <- 'integrated_snn_res.0.75'
Idents(seurat_integrated) <- res

######### VISUALISING MARKER GENES ------
DefaultAssay(seurat_integrated) <- "RNA"
markers <- read.xlsx("data/MANUSCRIPT/Drive Manuscript Supplementary Tables.xlsx", sheet = 2)[1:29,c(1:6)] 
row_order <- c("E", "M", "CYSC", "EC", "LC", "C", "G", "G, PS", "PS, SS", "PS, SS, ST", "PS, ST", "ST")
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
}

new_names <- sapply(mk_genes, new_name_func)

dps_all <- DotPlot(seurat_integrated, features = mk_genes, assay = "RNA")+coord_flip() + 
  scale_x_discrete(labels = new_names)
DP <- DimPlot(seurat_integrated, label = T, split.by = 'treatment')
FP <- FeaturePlot(seurat_integrated, features = 'nFeature_RNA', split.by = 'treatment')
FP2 <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.75, y = nFeature_RNA)) + geom_boxplot() + 
  facet_wrap(~treatment)

ggarrange(plotlist = list(dps_all, DP, FP, FP2), ncol = 4, nrow = 1, labels = c("A", "B", "C")) %>% 
  ggsave("plots/celltype_data_prep.pdf", plot = ., height = 10, width = 49)


samplects_counts <- seurat_integrated@meta.data %>% 
  group_by(sample, integrated_snn_res.0.75) %>% 
  summarise(n = n())
sample_counts <- seurat_integrated@meta.data %>% 
  group_by(sample) %>% summarise(n = n())

samplects_counts <- samplects_counts %>% 
  merge(sample_counts, by = "sample") %>% 
  mutate(freq = n.x/n.y) %>% 
  mutate(celltype = integrated_snn_res.0.75)

cellpopulation_compositions <- samplects_counts %>% dplyr::select(sample, freq, celltype) %>% 
  spread(., key = celltype, value = freq) %>% dplyr::select(-sample) %>% t() %>% 
  as.data.frame() %>% #add colnames as samples 
  `colnames<-`(sample_counts$sample) %>% #turn NAs into 0
  replace(is.na(.), 0)
cpcs <- 100 * cellpopulation_compositions/rowSums(cellpopulation_compositions)

which_ct_keep_function <- function(x){
  keep <- sum(x > 5) >= 3
  if (keep){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

keep <- apply(cpcs, 1, which_ct_keep_function)

SR_keep <- apply(cpcs[,1:4], 1, which_ct_keep_function)
ST_keep <- apply(cpcs[,5:8], 1, which_ct_keep_function)
keep <- SR_keep | ST_keep

#################### FINAL ASSIGNMENT -------------------------
seurat_integrated@meta.data$celltype <- "Keep"

#Limited marker expression, unknown celltypes
no_cell_cycle_marker_unknown_cells_0.75 <- c(1,7,15,23)
seurat_integrated@meta.data$celltype[seurat_integrated@meta.data$integrated_snn_res.0.75 %in%
                                       no_cell_cycle_marker_unknown_cells_0.75] <- "Unclear marker expression"
#Low representation
low_rep <- names(which(keep == F)) %>% as.numeric()
seurat_integrated@meta.data$celltype[seurat_integrated@meta.data$integrated_snn_res.0.75 %in% 
                                       low_rep] <- "Low representation"
seurat_integrated_ss <- subset(seurat_integrated, (celltype %in% c("Keep")))











###############################################################

#### #### #### #### #### ROUND 2 #### #### #### #### #### #### 

###############################################################

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
ggsave("plots/clustering_plots_ss.png", clustering_plots, height = 30, width = 30)

Idents(seurat_integrated_ss) <- 'integrated_snn_res.0.4'


######### VISUALISING MARKER GENES ------
DefaultAssay(seurat_integrated_ss) <- "RNA"

### Plotting key markers against numbered cell clusters for assignment ---
dps_all <- DotPlot(seurat_integrated_ss, features = mk_genes, assay = "RNA")+coord_flip() + 
  scale_x_discrete(labels = new_names)
DP <- DimPlot(seurat_integrated_ss, label = T, split.by = 'treatment')
FP <- FeaturePlot(seurat_integrated_ss, features = 'nFeature_RNA', split.by = 'treatment')
FP2 <- seurat_integrated_ss@meta.data %>% 
  ggplot(aes(x = integrated_snn_res.0.4, y = nFeature_RNA)) + geom_boxplot() + 
  facet_wrap(~treatment)
ggarrange(dps_all, DP, FP, FP2, ncol = 4) %>% 
  ggsave( "plots/subset_celltype_data_prep.pdf", ., width = 40, height = 10)
#--------------------------------------------------------------------------

#resolution 0.4 clusters
Muscle <- c(16)
Pre_meiotic_cyst <- c(2,5,6,8,15)
Post_meiotic_cyst <- c(10,3)
`GSC/Spermatogonia` <- c(0)
Primary_spermatocytes <- c(1,11)
Secondary_spermatocytes <- c(12,14)
Spermatids <- c(4,7,9,13,17)


#seurat_final <- CellSelector(DimPlot(seurat_final), seurat_final, "outliers1")
#seurat_final <- CellSelector(DimPlot(seurat_final), seurat_final, "outliers2")

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

all_figure_markers <- new_names
main_figure_markers <- new_names[]
save(seurat_final, main_figure_markers, all_figure_markers, file = 'data/RData/seurat_final.RData')







