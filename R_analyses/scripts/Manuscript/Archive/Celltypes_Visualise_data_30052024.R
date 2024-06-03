library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
load("data/RData/integrated_seurat_nf200_mtr0.20_gu0.RData")
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
pcs
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:pcs)
seurat_integrated <- RunPCA(seurat_integrated, dims = 1:pcs)


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
ggsave("plots/clustering_plots.png", clustering_plots, height = 30, width = 30)

res <- 'integrated_snn_res.0.4'
Idents(seurat_integrated) <- res

save(seurat_integrated, file = 'data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned.RData')

p1 <- seurat_integrated@meta.data %>% 
  dplyr::select(integrated_snn_res.0.4, treatment, sample) %>% 
  group_by(integrated_snn_res.0.4, treatment, sample) %>%
  summarise( # calculate number of cells in each cluster as proportion of total number of cells for SR and ST
    n = n()) %>% 
  group_by(treatment, sample) %>% 
  summarise(proportion = n/sum(n), 
            cluster = integrated_snn_res.0.4) %>% 
  ggplot(aes(x = cluster, y = proportion, fill = treatment)) + geom_boxplot()
  #ggplot(aes(x = cluster , y = n2, fill = treatment)) + geom_bar(position = "dodge", stat = "identity")
p2 <- DimPlot(seurat_integrated, group.by = "integrated_snn_res.0.4", label = T)
ggarrange(p1, p2, ncol = 1, nrow = 2)
######### VISUALISING MARKER GENES ------
DefaultAssay(seurat_integrated) <- "RNA"

markers_orths <- read.xlsx("data/MANUSCRIPT/DRIVE_MANUSCRIPT_SUPPLEMENTARY_TABLES.xlsx", sheet = 1)

marker_celltypes <- read.xlsx("data/markers/TDAL_CONSENSUS_MARKERS010324.xlsx", sheet = 2)[,c(1,2,11,13,14)]
markers <- markers_orths %>% 
  merge(marker_celltypes, by = c("Gene_name", "FBgn_name")) %>% 
  filter(!is.na(Confidence))

for (m in c("Consensus_Marker", "Medium_Group", "Broad_Group")){

  clusters <- str_split(markers[,m], ",") %>% unlist() %>% #remove space at start of string
    gsub("^ ", "", .) %>% unique()
  clusters <- clusters[!is.na(clusters)]
  #clusters <- markers$Medium_Group %>% unique()
  clusters <- markers[,m] %>% unique()
  
  clusters <- clusters[!is.na(clusters)]
  cluster_genes <- list()
  
  for (c in clusters){
    genes <- markers$T.dal_marker[markers[,m] == c]
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
  
  
  
  mk_genes <- markers$T.dal_marker %>% gsub(",", " ", .) %>% 
    str_split(" ") %>% unlist() %>% gsub("\\(|\\)", "", .) %>% 
    gsub("_", "-", .) %>% intersect(rownames(seurat_integrated@assays$RNA@counts))
  
  new_name_func <- function(gene){
    nn <- paste(markers_orths$Gene_name[grep(gene, gsub("_", "-", markers_orths$T.dal_marker ))], collapse = ", ")
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
#### HERE_ DOUBLE CHECK THE DATA TO MAKE SURE WHAT YOU DOWNLOADED IS THE SAME AND REDOWNLOAD TH SHEET IF NECESSARY
########## ASSIGN CELL TYPES HERE ----------------

GSC_G <- c(0,1)
Cyst <- c(3,5,10,11)
SR_specific <- c(4,6)
Early_Spermatid <- c(7)
Late_Spermatid <- c(8)
Primary_Spermatocyte <- c(13)
Secondary_Spermatocyte <- c(17)
Unknown <- c(2,9,12,14,15,16,18,19)

#make table with this information where column 1 is number and 2 is name of celltype
celltype_table <- data.frame(clus_no = c(0:19), celltype = NA)
celltype_table[celltype_table$clus_no %in% GSC_G, "celltype"] <- "GSC_G"
celltype_table[celltype_table$clus_no %in% Cyst, "celltype"] <- "Cyst"
celltype_table[celltype_table$clus_no %in% SR_specific, "celltype"] <- "SR_specific"
celltype_table[celltype_table$clus_no %in% Early_Spermatid, "celltype"] <- "Early_Spermatid"
celltype_table[celltype_table$clus_no %in% Late_Spermatid, "celltype"] <- "Late_Spermatid"
celltype_table[celltype_table$clus_no %in% Primary_Spermatocyte, "celltype"] <- "Primary_Spermatocyte"
celltype_table[celltype_table$clus_no %in% Secondary_Spermatocyte, "celltype"] <- "Secondary_Spermatocyte"
celltype_table[celltype_table$clus_no %in% Unknown, "celltype"] <- "Unknown"
celltype_table$celltype[celltype_table$celltype == "Unknown"] <- paste0("Unknown_", celltype_table[celltype_table$celltype == "Unknown", "clus_no"])

celltype_table$broad_celltype <- celltype_table$celltype
celltype_table$broad_celltype[celltype_table$clus_no %in% c(Early_Spermatid, Late_Spermatid)] <- "Spermatid"
celltype_table$broad_celltype[celltype_table$clus_no %in% c(Primary_Spermatocyte, Secondary_Spermatocyte)] <- "Spermatocyte"
celltype_table$germ_somatic <- "SOMA"
celltype_table$germ_somatic[celltype_table$clus_no %in% c(GSC_G, Early_Spermatid, Late_Spermatid, Primary_Spermatocyte, Secondary_Spermatocyte)] <- "GERM"
celltype_table$germ_somatic[celltype_table$clus_no %in% Unknown] <- "Unknown"

write.table(celltype_table, "data/celltype_table.txt", sep = "\t", row.names = F, col.names = T)

#----------------------------------------------

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


seurat_integrated@meta.data %>% 
  ggplot(aes(x = xprop, y = nFeature_SCT)) +
  geom_smooth(method = 'lm') + ggpubr::stat_cor() +
  geom_point()
DefaultAssay(seurat_integrated) <- "SCT"
seurat_integrated <- AverageExpression(seurat_integrated)

HtM <- DoHeatmap(seurat_integrated, features = unique(unlist(cluster_genes)), group.by = 'celltype')
ggsave("plots/tdal_markers_heatmap.pdf", height = 20, width = 20)

#################### FINAL ASSIGNMENT -------------------------

GSC_G <- c(0,1)
Cyst <- c(3,5,10,11)
SR_specific <- c(4,6)
Early_Spermatid <- c(7)
Late_Spermatid <- c(8)
Primary_Spermatocyte <- c(13)
Secondary_Spermatocyte <- c(17)
Unknown <- c(2,9,12,14,15,16,18,19)


seurat_integrated@metadata[[celltype]] <- NA
seurat_integrated@metadata[[celltype]][which(seurat_integrated@metadata$integrated_snn_res.0.4 %in% GSC_G)] <- "GSC_G"
seurat_integrated@metadata[[celltype]][which(seurat_integrated@metadata$integrated_snn_res.0.4 %in% Cyst)] <- "Cyst"
seurat_integrated@metadata[[celltype]][which(seurat_integrated@metadata$integrated_snn_res.0.4 %in% SR_specific)] <- "SR_specific"
seurat_integrated@metadata[[celltype]][which(seurat_integrated@metadata$integrated_snn_res.0.4 %in% Early_Spermatid)] <- "Early_Spermatid"
seurat_integrated@metadata[[celltype]][which(seurat_integrated@metadata$integrated_snn_res.0.4 %in% Late_Spermatid)] <- "Late_Spermatid"
seurat_integrated@metadata[[celltype]][which(seurat_integrated@metadata$integrated_snn_res.0.4 %in% Primary_Spermatocyte)] <- "Primary_Spermatocyte"
seurat_integrated@metadata[[celltype]][which(seurat_integrated@metadata$integrated_snn_res.0.4 %in% Secondary_Spermatocyte)] <- "Secondary_Spermatocyte"
seurat_integrated@metadata[[celltype]][which(seurat_integrated@metadata$integrated_snn_res.0.4 %in% Unknown)] <- "Unknown"

save(seurat_integrated, file = "data/RData//integrated_seurat_nf200_mtr0.20_gu0_cleaned_celltype.RData")

################################################################








########### BONUS MATERIAL --------------------

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
