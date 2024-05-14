###LIBRARIES ----
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)

files <- c("7-st1_stalkie_dros", "3-sr1_stalkie_dros")

#READING IN COUNT MATRICES and creating sample variables. 
sample_variables <- c()
for (file in files){
  print(file)
  var_name = str_split(file, "-", simplify = TRUE)[,2]
  seurat_data <- Read10X(data.dir = paste("indata/cellranger/unfiltered/", file, sep = ""))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 200, ### FILTERING CELLS WITH LESS THAN 200 GENES
                                   project = var_name) #Might be able to include min.cells = 3 for keeping a gene rather than dooing this later
  sample_variables <- c(sample_variables, var_name)
  assign(var_name, seurat_obj)
}

#Merging samples into single suerat object
merged_seurat <- merge(x = eval(parse(text = sample_variables[1])), 
                       c(y = eval(parse(text = sample_variables[2]))), 
                       add.cell.id = c("ST", "SR")) #ADDS IN SAMPLE TYPE/NAME TO CELL ID (ROWNAME IN METADATA)


merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA) # Calculating the number of features per UMI 

merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "mt")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100 #CALCULATING PROPORTION OF READS MAPPING TO MITOCHONDRIAL GENOME (x100 as originally a percentage)


metadata <- merged_seurat@meta.data # Seperately creating metadata dataframe to save metrics etc without risking affecting seurat object
metadata$cells <- rownames(metadata) #add cell IDs to metadat

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$sample <- metadata$seq_folder

merged_seurat@meta.data <- metadata #Save the more complete metadat to the seurat object


save(merged_seurat, file="outdata/seurat/merged_seurat.RData") #SAVING OUR OBJECT 

#FILTERING DATA
#Thresholds
#3/5/10 cells 
#100/200 features
#Remove with more than 5% transcriptome mitochodniral 
#
plot_data <- data.frame(seq_folder = NA, cells = NA, mean_nUMI = NA, 
                        mean_features = NA, kept_genes = NA, 
                        feature_filter = NA, cell_filter = NA, nUMI_filter = NA)

for (f in c(0, 100,200)){
  for (cells in c(0, 3,5,10)){
    for (umi in c(0, 100,250,350, 500)){
      filtered_seurat <- subset(x = merged_seurat, 
                                subset= (nUMI >= umi) & 
                                  (nGene >= f) & 
                                  (log10GenesPerUMI > 0) & # Can be dying cells or simple cell types such as bloody
                                  (mitoRatio < 0.05))
      counts <- GetAssayData(object = filtered_seurat, slot = "counts")
      count_func <- function(s, counts){
        sample_cells <- filtered_seurat@meta.data$cells[filtered_seurat@meta.data$seq_folder == s]
        ss <- counts[,sample_cells]
        kgs <- Matrix::rowSums(ss) >= cells #Filtering for genes that are expressed in at least 10 cells (thresholds = 1 read)
        return(names(kgs)[kgs])
      }
      
      keep_genes_per_sample <- sapply(unique(filtered_seurat@meta.data$seq_folder), count_func, counts = counts)
      keep_genes <- unique(c(unlist(keep_genes_per_sample)))
      filtered_counts <- counts[keep_genes, ]
      filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
      metadata_clean <- filtered_seurat@meta.data
      sum_data <-  metadata_clean %>% 
        group_by(seq_folder) %>% 
        summarise(cells = n(), 
                  mean_nUMI = mean(nUMI), 
                  mean_features = mean(nGene))
      sum_data$kept_genes <- nrow(filtered_counts)
      sum_data$feature_filter <- f
      sum_data$cell_filter <- cells
      sum_data$nUMI_filter <- umi
      plot_data <- rbind(plot_data, sum_data)
    }
  }
}

plot_data <- plot_data[-1,]
write.table(plot_data, "outdata/QC/paramater_checks.tsv", quote = FALSE, row.names = FALSE)
plots[[1]] <- ggplot(plot_data[plot_data$cell_filter == 0 & plot_data$feature_filter == 0,], aes(y = cells, x = as.factor(nUMI_filter), fill = seq_folder)) + geom_bar(stat = 'identity', position = 'dodge')
plots[[2]] <- ggplot(plot_data[plot_data$nUMI_filter == 0 & plot_data$feature_filter == 0,], aes(y = cells, x = as.factor(cell_filter), fill = seq_folder)) + geom_bar(stat = 'identity', position = 'dodge')
plots[[3]] <- ggplot(plot_data[plot_data$nUMI_filter == 0 & plot_data$cell_filter == 0,], aes(y = cells, x = as.factor(feature_filter), fill = seq_folder)) + geom_bar(stat = 'identity', position = 'dodge')
plots[[4]]<- ggplot(plot_data[plot_data$nUMI_filter == 0 & plot_data$feature_filter == 0,], aes(x = cell_filter, y = kept_genes, fill = seq_folder)) + geom_bar(stat = 'identity', position = 'dodge')
save_plots <- ggarrange(plotlist = plots)
ggsave('plots/comp_of_filters.pdf', save_plots, width = 10, height = 7)


#Filtering for lowly expressed genes across cells
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10 #Filtering for genes that are expressed in at least 10 cells (thresholds = 1 read)
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
metadata_clean <- filtered_seurat@meta.data


dim(metadata_clean[metadata_clean$seq_folder == "sr2_stalkie_dros",])


save(filtered_seurat, file="outdata/seurat/filtered_seurat.RData") #SAVING OUR OBJECT 



split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]])
}

integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000)

split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)


seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")



saveRDS(seurat_integrated, "results/integrated_seurat.rds")

seurat_integrated <- FindVariableFeatures(seurat_integrated, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
seurat_integrated <- ScaleData(seurat_integrated)

seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, 
                        dims = 1:40,
                        reduction = "pca")
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
UMAPa <- DimPlot(seurat_integrated, group.by = "sample", pt.size = 0.001) +  ggtitle("") +
  scale_color_manual(labels = c("st1", "st2", "st3", "st4"), values = c("#000000","#CC79A7", "#E69F00", "#56B4E9"))
UMAPb <- DimPlot(seurat_integrated, split.by = "sample")
UMAPc <- DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6, pt.size = 0.00001) +  theme(legend.position="none")


#togeth <- ggarrange(ggarrange(UMAPa, UMAPc, ncol = 2, labels = c("A", "B")), 
#          ggarrange(UMAPb, labels = "C") ,nrow = 2) 
togeth2 <- ggarrange(UMAPc, ncol = 1, labels = c("A"), font.label = c(size = 32))
ggsave("plots/SSE_UMAPAB.pdf", UMAPc, height = 8, width = 8) 


















del <- subset(x = merged_seurat, (nUMI >= 200) & 
                            (nGene >= 200) & 
                            (mitoRatio < 0.05) & 
                seq_folder == "st1_stalkie_dros")
dim(del@meta.data)
