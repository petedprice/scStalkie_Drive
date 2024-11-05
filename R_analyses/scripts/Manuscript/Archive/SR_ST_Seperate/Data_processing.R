library(tidyverse)
library(Seurat)
library(openxlsx)
library(stringr)
library(gridExtra)
library(ggpubr)
library(clustree)
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")

approach = "cluster"
if (approach == "re_integrate"){
  load("data/RData/param_checks//integrated_seurat_nf200_mtr0.20_gu0.RData")
  seurat_integrated$treatment <- "SR"
  seurat_integrated$treatment[grep("st", seurat_integrated$sample)] <- "ST"
  DefaultAssay(seurat_integrated) <- "RNA"
  seurat_integrated[['SCT']] <- NULL
  seurat_integrated[['integrated']] <- NULL
  
  
  split_seurat <- seurat_integrated %>% SplitObject(., split.by = "sample")
  
  #Prepossessing step necessary if SCT transformed
  split_seurat <- sapply(split_seurat, SCTransform, vars.to.regress =
                           c("mitoRatio","nUMI","S.Score","G2M.Score"))
  
  samples <- list(c("sr1", "sr2", "sr3", "sr5"), c("st1", "st2", "st3", "st5"))
  names(samples) <- c("SR", "ST")
  
  SR_features <- SelectIntegrationFeatures(object.list = split_seurat[samples[[1]]], nfeatures = 3000)
  ST_features <- SelectIntegrationFeatures(object.list = split_seurat[samples[[2]]], nfeatures = 3000)
  
  SR_split_seurat <- PrepSCTIntegration(object.list = split_seurat[samples[[1]]],
                                     anchor.features = SR_features)
  ST_split_seurat <- PrepSCTIntegration(object.list = split_seurat[samples[[2]]],
                                     anchor.features = ST_features)
  
  SR_anchors <- FindIntegrationAnchors(object.list = SR_split_seurat,
                                    anchor.features = SR_features,
                                    normalization.method = "SCT")
  ST_anchors <- FindIntegrationAnchors(object.list = ST_split_seurat,
                                    anchor.features = ST_features,
                                    normalization.method = "SCT")
  
  SR_seurat_integrated <- IntegrateData(anchorset = SR_anchors,
                                        normalization.method = "SCT")
  
  ST_seurat_integrated <- IntegrateData(anchorset = ST_anchors,
                                        normalization.method = "SCT")
  
  seur_objs <- list(SR_seurat_integrated, ST_seurat_integrated)
  names(seur_objs) <- c("SR", "ST")
  
  ######### Determining number of PCAs to use for clustering/UMAP etc.
  DefaultAssay(seur_objs[[1]]) <- "integrated"
  DefaultAssay(seur_objs[[2]]) <- "integrated"
  for (i in 1:length(seur_objs)){
    seur_objs[[i]] <- RunPCA(seur_objs[[i]], dims = 1:40)
    EP <- ElbowPlot(seur_objs[[i]], ndims = 40)
    
    pct <- seur_objs[[i]][["pca"]]@stdev / sum(seur_objs[[i]][["pca"]]@stdev) * 100
    
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    
    # last point where change of % of variation is more than 0.1%.
    pcs <- min(co1, co2)
    seur_objs[[i]] <- RunUMAP(seur_objs[[i]], dims = 1:pcs)
    seur_objs[[i]] <- RunPCA(seur_objs[[i]], dims = 1:pcs)
    seur_objs[[i]] <- RunTSNE(seur_objs[[i]], dims = 1:pcs)
    
    ######### Determining resolution to use --------
    seur_objs[[i]] <- FindNeighbors(object=seur_objs[[i]], dims=1:pcs)
    seur_objs[[i]] <- FindClusters(object=seur_objs[[i]], resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                                                               0.75, 1, 1.25, 1.5, 1.75, 2, 
                                                                               2.5, 3))
    clusttree <- clustree(seur_objs[[i]])
    CT <- plot(clusttree)
    CT_noclusters <- clusttree$data %>% 
      group_by(integrated_snn_res.) %>% 
      summarise(no_clusters = n_distinct(cluster)) %>% 
      dplyr::rename(resolution = `integrated_snn_res.`) %>%
      ggplot(aes(x = resolution, y = no_clusters)) +
      geom_point()
    
    clustering_plots <- ggarrange(EP, CT, CT_noclusters, ncol = 1, nrow = 3)
    ggsave(paste0("plots/", names(seur_objs)[i], "_clustering_plots.png"), clustering_plots, height = 30, width = 30)
    res <- 'integrated_snn_res.0.4'
    Idents(seur_objs[[i]]) <- res
  }
  save(seur_objs, file = "data/RData/seur_objs_integrated.RData")
  
} else if (approach = "marker_checks"){
  
  #LOAD DATA 
  load("data/RData/seur_objs_integrated.RData")
  
  ######### VISUALISING MARKER GENES ------
  markers <- read.xlsx("data/MANUSCRIPT/DRIVE_MANUSCRIPT_SUPPLEMENTARY_TABLES.xlsx", sheet = 1)[,c(2,3,4,5,6,7)] %>% 
    filter(`DOI.(Source)` != "https://doi.org/10.1038/s41467-021-20897-y")
  for (i in 1:length(seur_objs)){
    
    DefaultAssay(seur_objs[[i]]) <- "RNA"
    
    m <- "Celltype.(Specific)"
    clusters <- markers[,m] %>% unique()
    clusters <- clusters[!is.na(clusters)]
    cluster_genes <- list()
    for (c in clusters){
      genes <- markers$T.dal.Ortholog[markers[,m] == c]
      genes2 <- genes %>% gsub(",", " ", .) %>% 
        str_split(" ") %>% unlist() %>% gsub("\\(|\\)", "", .) %>% 
        gsub("_", "-", .) %>% intersect(rownames(seur_objs[[i]]@assays$RNA@counts))
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
      gsub("_", "-", .) %>% intersect(rownames(seur_objs[[i]]@assays$RNA@counts))
    
    new_name_func <- function(gene){
      nn <- paste(markers$`Drosophila.Marker.(Gene.Name)`[grep(gene, gsub("_", "-", markers$T.dal.Ortholog ))], collapse = ", ")
      nn2 <- paste0(nn , " (", gene, ")")
    }
    
    new_names <- sapply(mk_genes, new_name_func)
    dps_list <- lapply(cluster_genes, function(x) DotPlot(seur_objs[[i]], features = x, assay = "RNA")+coord_flip() +   
                         scale_x_discrete(labels = new_names))
    dps_all <- DotPlot(seur_objs[[i]], features = unique(unlist(one_cluster_genes)), assay = "RNA")+coord_flip() + 
      scale_x_discrete(labels = new_names)
    dps <- ggarrange(plotlist = dps_list, ncol = 2, nrow = ceiling(length(dps_list)/2), labels = names(cluster_genes))
    
    #arrange them again but with names as titles
    ggsave(paste0("plots/", names(seur_objs)[i], "_tdal_markers_dotplot_clus.png"), dps, height = 30, width = 30)
    ggsave(paste0("plots/", names(seur_objs)[i], "_tdal_markers_dotplot_all_clus.pdf"), dps_all, height = 20, width = 14)
    dmp <- DimPlot(seur_objs[[i]], label = T)
    ggsave(paste0("plots/", names(seur_objs)[i], "_tdal_markers_dimplot_clus.pdf"), dmp, height = 10, width = 10)
    
    
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
    dps_all_figure <- DotPlot(seur_objs[[i]], features = names(main_figure_markers), assay = "RNA")+coord_flip() +
      scale_x_discrete(labels = as.vector((main_figure_markers)))
    DP <- DimPlot(seur_objs[[i]], label = T)
    FP <- FeaturePlot(seur_objs[[i]], features = 'nFeature_RNA')
    FP2 <- seur_objs[[i]]@meta.data %>% 
      ggplot(aes(x = integrated_snn_res.0.4, y = nFeature_RNA)) + geom_boxplot()
    ggarrange(plotlist = list(dps_all_figure, DP, FP, FP2), ncol = 4, nrow = 1, labels = c("A", "B", "C")) %>% 
      ggsave(paste0("plots/", names(seur_objs)[i], "_main_figure.pdf"), plot = ., height = 10, width = 40)
  }
} else if (approach = "subsetting"){
  #################### FINAL ASSIGNMENT -------------------------
  #For No max feature filter 
  ST_Unknown <-c(5,8,13,14,15,17)
  SR_Unknown <- c(2,6,12,15)
  
  seur_objs$ST@meta.data$celltype <- "Keep"
  seur_objs$SR@meta.data$celltype <- "Keep"
  
  
  seur_objs$ST@meta.data$celltype[which(seur_objs$ST@meta.data$integrated_snn_res.0.4 %in% ST_Unknown)] <- "Unknown"
  seur_objs$SR@meta.data$celltype[which(seur_objs$SR@meta.data$integrated_snn_res.0.4 %in% SR_Unknown)] <- "Unknown"
  
  seur_objs_ss <- seur_objs
  seur_objs_ss$ST <- subset(seur_objs_ss$ST, celltype == "Keep")
  seur_objs_ss$SR <- subset(seur_objs_ss$SR, celltype == "Keep")
  
  save(seur_objs_ss, file = "data/RData/re_integrated_seurat_objs_ss.RData")
###############################################################
} else if (approach == "reclustering"){
  load("data/RData/re_integrated_seurat_objs_ss.RData")
  for (i in 1:length(seur_objs_ss)){
    DefaultAssay(seur_objs_ss[[i]]) <- "integrated"
    seur_objs_ss[[i]] <- RunPCA(seur_objs_ss[[i]], dims = 1:40)
    EP <- ElbowPlot(seur_objs_ss[[i]], ndims = 40)
    
    pct <- seur_objs_ss[[i]][["pca"]]@stdev / sum(seur_objs_ss[[i]][["pca"]]@stdev) * 100
    
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    
    # last point where change of % of variation is more than 0.1%.
    pcs <- min(co1, co2)
    seur_objs_ss[[i]] <- RunUMAP(seur_objs_ss[[i]], dims = 1:pcs)
    seur_objs_ss[[i]] <- RunPCA(seur_objs_ss[[i]], dims = 1:pcs)
    seur_objs_ss[[i]] <- RunTSNE(seur_objs_ss[[i]], dims = 1:pcs)
    
    ######### Determining resolution to use --------
    seur_objs_ss[[i]] <- FindNeighbors(object=seur_objs_ss[[i]], dims=1:pcs)
    seur_objs_ss[[i]] <- FindClusters(object=seur_objs_ss[[i]], resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                                                         0.75, 1, 1.25, 1.5, 1.75, 2, 
                                                                         2.5, 3))
    clusttree <- clustree(seur_objs_ss[[i]])
    CT <- plot(clusttree)
    CT_noclusters <- clusttree$data %>% 
      group_by(integrated_snn_res.) %>% 
      summarise(no_clusters = n_distinct(cluster)) %>% 
      dplyr::rename(resolution = `integrated_snn_res.`) %>%
      ggplot(aes(x = resolution, y = no_clusters)) +
      geom_point()
    
    clustering_plots <- ggarrange(EP, CT, CT_noclusters, ncol = 1, nrow = 3)
    ggsave(paste0("plots/", names(seur_objs_ss)[i], "_clustering_plots_ss.png"), clustering_plots, height = 30, width = 30)
    res <- 'integrated_snn_res.1'
    Idents(seur_objs_ss[[i]]) <- res
  }
  
} else if (approach == "final"){
  ######### VISUALISING MARKER GENES ------
  markers <- read.xlsx("data/MANUSCRIPT/DRIVE_MANUSCRIPT_SUPPLEMENTARY_TABLES.xlsx", sheet = 1)[,c(2,3,4,5,6,7)] %>% 
    filter(`DOI.(Source)` != "https://doi.org/10.1038/s41467-021-20897-y")
  for (i in 1:length(seur_objs_ss)){
    
    DefaultAssay(seur_objs_ss[[i]]) <- "RNA"
    
    m <- "Celltype.(Specific)"
    clusters <- markers[,m] %>% unique()
    clusters <- clusters[!is.na(clusters)]
    cluster_genes <- list()
    for (c in clusters){
      genes <- markers$T.dal.Ortholog[markers[,m] == c]
      genes2 <- genes %>% gsub(",", " ", .) %>% 
        str_split(" ") %>% unlist() %>% gsub("\\(|\\)", "", .) %>% 
        gsub("_", "-", .) %>% intersect(rownames(seur_objs_ss[[i]]@assays$RNA@counts))
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
      gsub("_", "-", .) %>% intersect(rownames(seur_objs_ss[[i]]@assays$RNA@counts))
    
    new_name_func <- function(gene){
      nn <- paste(markers$`Drosophila.Marker.(Gene.Name)`[grep(gene, gsub("_", "-", markers$T.dal.Ortholog ))], collapse = ", ")
      nn2 <- paste0(nn , " (", gene, ")")
    }
    
    new_names <- sapply(mk_genes, new_name_func)
    dps_list <- lapply(cluster_genes, function(x) DotPlot(seur_objs_ss[[i]], features = x, assay = "RNA")+coord_flip() +   
                         scale_x_discrete(labels = new_names))
    dps_all <- DotPlot(seur_objs_ss[[i]], features = unique(unlist(one_cluster_genes)), assay = "RNA")+coord_flip() + 
      scale_x_discrete(labels = new_names)
    dps <- ggarrange(plotlist = dps_list, ncol = 2, nrow = ceiling(length(dps_list)/2), labels = names(cluster_genes))
    
    #arrange them again but with names as titles
    ggsave(paste0("plots/", names(seur_objs_ss)[i], "_tdal_markers_dotplot_clus_ss.png"), dps, height = 30, width = 30)
    ggsave(paste0("plots/", names(seur_objs_ss)[i], "_tdal_markers_dotplot_all_clus_ss.pdf"), dps_all, height = 20, width = 14)
    dmp <- DimPlot(seur_objs_ss[[1]], label = T)
    ggsave(paste0("plots/", names(seur_objs_ss)[i], "_tdal_markers_dimplot_clus_ss.pdf"), dmp, height = 10, width = 10)
    
    
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
    dps_all_figure <- DotPlot(seur_objs_ss[[i]], features = names(main_figure_markers), assay = "RNA")+coord_flip() +
      scale_x_discrete(labels = as.vector((main_figure_markers)))
    DP <- DimPlot(seur_objs_ss[[i]], label = T)
    FP <- FeaturePlot(seur_objs_ss[[i]], features = 'nFeature_RNA')
    FP2 <- seur_objs_ss[[i]]@meta.data %>% 
      ggplot(aes(x = integrated_snn_res.0.4, y = nFeature_RNA)) + geom_boxplot()
    ggarrange(plotlist = list(dps_all_figure, DP, FP, FP2), ncol = 4, nrow = 1, labels = c("A", "B", "C")) %>% 
      ggsave(paste0("plots/", names(seur_objs_ss)[i], "_main_figure_ss.pdf"), plot = ., height = 10, width = 40)
    
  }
  
  seur_objs_ss$SR$celltype[seur_objs_ss$SR$integrated_snn_res.0.4 %in% c(13)] <- "Muscle"
  seur_objs_ss$SR$celltype[seur_objs_ss$SR$integrated_snn_res.0.4 %in% c(5,9,11,6)] <- "Pre-meiotic cyst"
  seur_objs_ss$SR$celltype[seur_objs_ss$SR$integrated_snn_res.0.4 %in% c(4,7)] <- "Post-meiotic cyst"
  seur_objs_ss$SR$celltype[seur_objs_ss$SR$integrated_snn_res.0.4 %in% c(2,8)] <- "GSC/Spermatogonia"
  seur_objs_ss$SR$celltype[seur_objs_ss$SR$integrated_snn_res.0.4 %in% c(12)] <- "Primary Spermatocytes"
  seur_objs_ss$SR$celltype[seur_objs_ss$SR$integrated_snn_res.0.4 %in% c(10)] <- "Secondary Spermatocytes"
  seur_objs_ss$SR$celltype[seur_objs_ss$SR$integrated_snn_res.0.4 %in% c(0,1,3)] <- "Spermatids"
  
  #now for ST
  seur_objs_ss$ST$celltype[seur_objs_ss$ST$integrated_snn_res.0.4 %in% c(11)] <- "Muscle"
  seur_objs_ss$ST$celltype[seur_objs_ss$ST$integrated_snn_res.0.4 %in% c(12,5,10,9)] <- "Pre-meiotic cyst"
  seur_objs_ss$ST$celltype[seur_objs_ss$ST$integrated_snn_res.0.4 %in% c(6)] <- "Post-meiotic cyst"
  seur_objs_ss$ST$celltype[seur_objs_ss$ST$integrated_snn_res.0.4 %in% c(2,4,1)] <- "GSC/Spermatogonia"
  seur_objs_ss$ST$celltype[seur_objs_ss$ST$integrated_snn_res.0.4 %in% c(13,8)] <- "Primary Spermatocytes"
  seur_objs_ss$ST$celltype[seur_objs_ss$ST$integrated_snn_res.0.4 %in% c(3)] <- "Secondary Spermatocytes"
  seur_objs_ss$ST$celltype[seur_objs_ss$ST$integrated_snn_res.0.4 %in% c(0,7)] <- "Spermatids"
  
  Idents(seur_objs_ss$SR) <- "celltype"
  Idents(seur_objs_ss$ST) <- "celltype"
  
  levels(seur_objs_ss$SR) <- c("Muscle", "Pre-meiotic cyst", "Post-meiotic cyst", 
                            "GSC/Spermatogonia", "Primary Spermatocytes", 
                            "Secondary Spermatocytes", "Spermatids")
    
  levels(seur_objs_ss$ST) <- c("Muscle", "Pre-meiotic cyst", "Post-meiotic cyst", 
                            "GSC/Spermatogonia", "Primary Spermatocytes", 
                            "Secondary Spermatocytes", "Spermatids")
  
  save(seur_objs_ss, file = "data/RData/seur_objs_ss_integrated_celltypes.RData")
}








