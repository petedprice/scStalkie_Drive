library(tidyverse)
library(Seurat)
library(ggpubr)
library(ggExtra)
library(patchwork)

seur_objs <- list.files("data/RData/param_checks/", full.names = TRUE, pattern = "seurat")

check_values_function <- function(obj, plot = F){
  outdata <- list()
  load(obj)
  seurat_integrated$treatment <- "ST"
  seurat_integrated$treatment[grepl("^sr", seurat_integrated$sample)] <- "SR"
  
  seurat_integrated$celltype <- seurat_integrated$integrated_snn_res.0.1
  
  sample_counts <- table(seurat_integrated$sample)
  sample_counts <- rbind(names(sample_counts), sample_counts, rep("ST"))
  sample_counts[3,grep("sr", colnames(sample_counts))] <- "SR"
  sc2 <-  cbind(t(sample_counts), seurat_integrated@meta.data[1,c("featurefilter", "mitofilter", "complexityfilter")])
  colnames(sc2)[1:3] <- c('sample', 'freq', "treatment")
  
  outdata$metadata <- seurat_integrated@meta.data
  outdata$samples <- sc2
  if (plot == T){
    title = paste("feature_filt >", sc2$featurefilter[1], ",mito_filt <", 
                  sc2$mitofilter[1], ",complx_filt >", sc2$complexityfilter[1])
    DefaultAssay(seurat_integrated) <- "RNA"
    FeaturePlots <- list()
    FeaturePlots$mitoRatio <- FeaturePlot(seurat_integrated, features = "mitoRatio", split.by = 'treatment') & 
      theme(legend.position = c(0.1,0.2))
    FeaturePlots$log10GenesPerUMI <- FeaturePlot(seurat_integrated, features = "log10GenesPerUMI", split.by = 'treatment')& 
      theme(legend.position = c(0.1,0.2))
    FeaturePlots$nFeature_RNA <- FeaturePlot(seurat_integrated, features = "nFeature_RNA", split.by = 'treatment')& 
      theme(legend.position = c(0.1,0.2))
    FeaturePlots$nCount_RNA <- FeaturePlot(seurat_integrated, features = "nCount_RNA", split.by = 'treatment')& 
      theme(legend.position = c(0.1,0.2))
    FeaturePlots$DimPlot <- DimPlot(seurat_integrated, group.by = 'celltype', split.by = 'treatment', label = T) & 
      theme(legend.position = c(0.1,0.2))
    FeaturePlots <- lapply(FeaturePlots, function(x)(return(x + plot_annotation(title = title))))
    outdata$FeaturePlots <- FeaturePlots
  }
  return(outdata)
}

allplots <- lapply(seur_objs, check_values_function, plot = T)
names(allplots) <- seur_objs
cell_numbers <- allplots %>% 
  lapply(., function(x)(return(x$samples))) %>% 
  bind_rows() %>% 
  mutate(freq = as.numeric(freq))

metadata <- allplots %>% 
  lapply(., function(x)(return(x$metadata))) %>% 
  bind_rows()

#see relationship between mito and complexity. 
#- A plot for the unfiltered data, plot mitoRatio against the log10genes/umi value
#- Here include the FeautrePlot of unfiltered mitoRatio and unfiltered log10genes/umi to see how it’s distributed
mito_against_gu_plots <- list()
mito_against_gu_plots[[1]] <- metadata %>% filter(featurefilter == 100, mitofilter == 1, complexityfilter == 0) %>% 
  ggplot(aes(x = mitoRatio, y = log10GenesPerUMI, colour = treatment)) + geom_point(alpha = 0.1) #add density plots onto x and y axis
mito_against_gu_plots[[1]] <- ggMarginal(mito_against_gu_plots[[1]], 
                                         groupColour = T, groupFill = T)
mito_against_gu_plots[[2]] <- metadata %>% filter(featurefilter == 100, mitofilter == 1, complexityfilter == 0) %>%
  ggplot(aes(x = mitoRatio, y = nFeature_RNA, colour = treatment)) + geom_point(alpha = 0.1) #add density plots onto x and y axis
mito_against_gu_plots[[2]] <- ggMarginal(mito_against_gu_plots[[2]], 
                                         groupColour = T, groupFill = T)
mito_against_gu_plots[[3]] <- metadata %>% filter(featurefilter == 100, mitofilter == 1, complexityfilter == 0) %>%
  ggplot(aes(x = mitoRatio, y = log(nCount_RNA), colour = treatment)) + geom_point(alpha = 0.1) #add density plots onto x and y axis
mito_against_gu_plots[[3]] <- ggMarginal(mito_against_gu_plots[[3]], 
                                         groupColour = T, groupFill = T)

mito_against_gu_plots[[4]] <- allplots$
  `data/RData/param_checks//integrated_seurat_nf100_mtr1_gu0.RData`$
  FeaturePlots$mitoRatio
mito_against_gu_plots[[5]] <- allplots$
  `data/RData/param_checks//integrated_seurat_nf100_mtr1_gu0.RData`$
  FeaturePlots$log10GenesPerUMI


mito_against_gu_plots[[6]] <- metadata %>% filter(featurefilter == 100, mitofilter == 1, complexityfilter == 0) %>%
  ggplot(aes(y = mitoRatio, fill = treatment,  x = celltype)) + geom_boxplot() #add density plots onto x and y axis
mito_against_gu_plots[[7]] <- metadata %>% filter(featurefilter == 100, mitofilter == 1, complexityfilter == 0) %>%
  ggplot(aes(y = log10GenesPerUMI, fill = treatment,  x = celltype)) + geom_boxplot() #add density plots onto x and y axis

mito_against_gu_plots[[8]] <- allplots$
  `data/RData/param_checks//integrated_seurat_nf100_mtr1_gu0.RData`$
  FeaturePlots$DimPlot


magp_arranged <- ggarrange(plotlist = mito_against_gu_plots)
ggsave("plots/magp_arranged.pdf", magp_arranged, width = 20, height = 20)

#Number of cells before and after filtering for mito expression 
#- So just do a table of boxplot for each treatment with mito filtering on the X and number of cells on Y, then have the complexity and nfeatures as X and Y 


cellnumbers_plots <- cell_numbers %>% 
  ggplot(aes(x = as.factor(mitofilter), y = freq, fill = treatment)) + geom_boxplot() + 
  facet_wrap(~complexityfilter+featurefilter)

cell_numbers_table <- cell_numbers %>% 
  group_by(treatment, mitofilter, complexityfilter, featurefilter) %>%
  mutate(freq = as.numeric(freq)) %>%
  summarise(freq_total = sum(freq), 
            mean_freq = mean(freq),
            sd_freq = sd(freq)) %>% 
  ggtexttable() 

ggarrange(plotlist = list(cellnumbers_plots, cell_numbers_table), ncol = 2, widths = c(3,3))
ggsave("plots/cell_numbers.pdf", width = 14, height = 10)

