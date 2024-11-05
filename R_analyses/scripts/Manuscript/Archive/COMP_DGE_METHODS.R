################# load libraries and functions and data prep ###################
library(dplyr)
library(Seurat)
library(openxlsx)
library(ggpubr)
library(scuttle)
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(ggrepel)
library(lme4)
library(GenomicFeatures)
library(ggrepel)
library(lme4)


load("data/RData/seurat_final.RData")
seurat_final$celltype
DefaultAssay(seurat_final) <- "RNA"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

Xgenes <- gsub("_", "-", filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME) %>% 
  intersect(rownames(seurat_final))
Agenes <- gsub("_", "-", filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME) %>% 
  intersect(rownames(seurat_final))
sce <- seurat_final %>% as.SingleCellExperiment(assay = "RNA")
dim(sce[rowSums(counts(sce)>1)>=10,])

################################################################################



####################### PLOT PCAS and DGE FUNCTION #############################
make_pca_function <- function(x, ortholog_table){
  mds_data <- x$logcpm %>% 
    plotMDS()
  
  PCA <- data.frame(x = mds_data$x, 
                    y = mds_data$y, 
                    names = x$res$samples$seq_folder,
                    treatment = x$res$samples$treatment) %>% 
    ggplot(aes(x = x, y = y, fill = names, label = names, colour = treatment))  +
    geom_text(hjust="inward", vjust="inward") +
    labs(x = "PC2", y = "PC1", title = x$res$samples$celltype[1]) 
  
  res <- x$toptags$table[is.na(x$toptags$table$logFC) == FALSE,] %>%
    as.data.frame()
  res$geneID <- rownames(res)
  
  res <- merge(res, ortholog_table, by.x = 'geneID', by.y = 'REF_GENE_NAME')   
  top_p <- res$FDR[order(res$FDR)][5]
  top_p2 <- res$FDR[order(res$FDR)][5]
  res$celltype = x$res$samples$celltype[1]
  Volcano <- res %>% 
    ggplot(aes(x = logFC, y = -log10(FDR), label = consensus_gene,
               colour = (abs(logFC) > 1 &  FDR < 0.05))) + 
    geom_point() + geom_text_repel(data=subset(res, FDR < top_p & abs(logFC) > 1),
                                   aes(logFC,-log10(FDR),label=consensus_gene), 
                                   hjust="inward", vjust="inward") +
    labs(y = "-log10 Pvalue (FDR adjusted)", x = "log2 fold-change", 
         title = paste(x$res$samples$celltype[1], "(positive is ST bias)"))+ 
    theme(legend.position = "none")  + 
    geom_hline(yintercept = -log10(0.05), color = 'grey', linetype = 'dashed') + 
    geom_vline(xintercept = c(-1,1), color = 'grey', linetype = 'dashed')
  ST_bias_genes <- filter(res, logFC > 1, FDR < 0.05)
  SR_bias_genes <- filter(res, logFC < -1, FDR < 0.05)
  output <- list(PCA,Volcano, ST_bias_genes, SR_bias_genes, res)
  names(output) <- c("PCA", "Volcano", "ST_bias", "SR_bias", "all")
  return(output)
  
}

################################################################################


######################## PSEUDOBULK INDIVIDUAL CELL TYPES ######################
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]
y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))
design <- model.matrix(~0 + treatment + celltype, data = agg_cell_types@colData)


keep <- filterByExpr(y, group = paste0(y$samples$celltype, y$samples$treatment), 
                     large.n = 3, min.prop = 0.49, min.count = 5)

y <- y[keep,]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
cpm <- edgeR::cpm(y, log=T) %>% as.data.frame()
dim(cpm)
mds_data <-cpm %>% 
  plotMDS()
plot <- data.frame(x = mds_data$x, 
                   y = mds_data$y, 
                   names = y$samples$seq_folder,
                   treatment = y$samples$treatment, 
                   cell_type = y$samples$celltype) %>% 
  ggplot(aes(x = x, y = y, label = cell_type, colour = treatment))  +
  geom_text(hjust="inward", vjust="inward") +
  labs(x = "PC2", y = "PC1", title = "ind cell types")

pdf("plots/PCA_pseudobulk_ind_cell_types.pdf")
plot 
dev.off()


################################################################################


## NEW PBDGE FUNCTION 


PBDGE <- function(label, sce){
  dge_data <- sce[,sce$celltype == label]
  y <- DGEList(counts(dge_data), samples=colData(dge_data), group = dge_data$treatment)
  design <- model.matrix(~treatment, y$samples)
  y <- normLibSizes(y)
  y <- estimateDisp(y, design)
  cpm <- cpm(y, log = F)
  
  keep1 <- cpm[,y$samples$treatment== "ST"] >= 2
  keep1 <- rowSums(keep1) > (ncol(keep1)/2)
  keep2 <- cpm[,y$samples$treatment== "SR"] >= 2
  keep2 <- rowSums(keep2) > (ncol(keep2)/2)
  keep <- keep1 | keep2
  
  #keep <- filterByExpr(y, group = y$sample$treatment, 
  #                     large.n = 3, min.prop = 0.5, min.count = 5)
  y <- y[keep,]
  
  STSR <- exactTest(y, pair = c("ST", "SR")) 
  fit <- glmQLFit(y, design, robust = T)
  res <- glmQLFTest(fit, coef=ncol(design))
  out <- list()
  out$res <- res
  out$logcpm <- cpm(y, log = T, prior.count = .0001)
  out$cpm <- cpm(y, log = F)
  out$toptags <- topTags(res, n = nrow(res))
  out$exacttest <- topTags(STSR, n = nrow(STSR))
  return(out)
}

############################## CELLTYPE DGE ------------------------------------

de.results_act <- pseudoBulkDGE(agg_cell_types, 
                                label=agg_cell_types$celltype,
                                condition = agg_cell_types$treatment,
                                design=~treatment,
                                coef = "treatmentST"
)


de.results_act <- sapply(unique(agg_cell_types$celltype), PBDGE, sce = agg_cell_types, 
                         USE.NAMES = T, simplify = F)
names(de.results_act) <- unique(agg_cell_types$celltype)


cell_type_DEG_output <- lapply(de.results_act, make_pca_function, ortholog_table = ortholog_table)


PCAs <- lapply(cell_type_DEG_output, function(x)(return(x$PCA)))
Volcanos <- lapply(cell_type_DEG_output, function(x)(return(x$Volcano)))

pdf("plots/PCA_compiled.pdf", width = 14, height = 10)
ggarrange(plotlist =  PCAs, common.legend = T)
dev.off()

pdf("plots/Volcanos_compiled.pdf", width = 20, height = 15)
ggarrange(plotlist = Volcanos, common.legend = T)
dev.off()

################################################################################




########################## DOSAGE COMPENSATION PLOTS ---------------------------

raw_counts <- y$counts %>% as.data.frame()
cpm$genes <- rownames(cpm)
raw_counts$genes <- rownames(raw_counts)

cols <- c("treatment", "celltype", "sample")
tidy_cpm <- cpm %>% as.data.frame %>% 
  pivot_longer(cols = colnames(cpm)[1:(length(colnames(cpm))-1)], 
               names_to = "sample", values_to = "logcpm") %>% 
  merge(y$samples[,cols], by.x = "sample", by.y = 0) %>% 
  rename(edgeRsample = sample,  sample = sample.y) %>% 
  mutate(Chr = ifelse(genes %in% Xgenes, "X", "A")) %>% 
  mutate(Chr = ifelse(!(genes %in% c(Xgenes, Agenes)), "Mito", Chr))

write.csv(tidy_cpm, "data/MANUSCRIPT/all_cpm_values.csv", row.names = F)
tidy_raw_counts <- raw_counts %>% as.data.frame %>% 
  pivot_longer(cols = colnames(raw_counts)[1:(length(colnames(raw_counts))-1)], 
               names_to = "sample", values_to = "raw_counts") %>% 
  merge(y$samples[,cols], by.x = "sample", by.y = 0) %>% 
  rename(edgeRsample = sample,  sample = sample.y) %>% 
  mutate(Chr = ifelse(genes %in% Xgenes, "X", "A")) %>% 
  mutate(Chr = ifelse(!(genes %in% c(Xgenes, Agenes)), "Mito", Chr))

################################################################################



############ DC 2 ------------

get_cpm_values_func <- function(x){
  md <- x$res
  ST <- rownames(filter(md$samples, treatment == "ST"))
  SR <- rownames(filter(md$samples, treatment == "SR"))
  cpm <- x$logcpm %>% 
    as.data.frame() %>%
    rownames_to_column("genes") %>% 
    pivot_longer(cols = colnames(.)[2:ncol(.)], names_to = "sample", values_to = "logcpm") %>% 
    mutate(Chr = ifelse(genes %in% Xgenes, "X", 
                        ifelse(genes %in% Agenes, "A", "Mito"))) %>% 
    mutate(treatment = ifelse(sample %in% ST, "ST", "SR")) %>% 
    mutate(celltype = md$samples$celltype[1]) %>% 
    rename(edgeRsample = sample) %>% 
    mutate(sample = md$samples[edgeRsample,]$sample)
  return(cpm)
}  

tidy_cpm2 <- lapply(de.results_act, get_cpm_values_func) %>% bind_rows()

write.csv(tidy_cpm2, "data/MANUSCRIPT/filtered_cpm_values.csv", row.names = F)
################################







######################## CHISQUARED ENRICHMENT FIGURE ##########################

dif_exp_data <- lapply(cell_type_DEG_output, function(x)(return(x$all))) %>% 
  bind_rows()

dif_exp_data <-  dif_exp_data %>% mutate(Significant = ifelse(dif_exp_data$logFC > 1 & dif_exp_data$FDR < 0.05, "ST-biased",
                                                              ifelse(dif_exp_data$logFC < -1 & dif_exp_data$FDR < 0.05, "SR-biased", "Unbiased")))


dif_exp_data_table <- dif_exp_data %>% 
  mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  group_by(celltype, Significant, chr) %>% 
  mutate(celltype = gsub("\n", "", celltype)) %>% 
  summarise(count = n()) %>% 
  reshape2::dcast(celltype + chr ~ Significant, value.var = "count") %>% 
  rename(`Cell type` = celltype, Chromosome = chr) %>% 
  .[c(3,4,7,8,5,6,1,2,9:14),]
dif_exp_data_table[is.na(dif_exp_data_table)] <- 0


### CHISQURED ---                                                                  
csq_chr_func <- function(x){
  data <- x$all
  data$sig <- 'Unbiased'
  data$sig[data$PValue < 0.05 & data$logFC > 1] <- 'ST-biased'
  data$sig[data$PValue < 0.05 & data$logFC < -1] <- 'SR-biased'
  chsq_table <- data %>%
    group_by(sig, chr) %>% 
    summarise(n = n()) %>% #trnasformn variables into rows and columns
    spread(sig, n)
  chsq_results <- chisq.test(chsq_table[,c(2:4)])
  rsds <- chsq_results$residuals %>% as.data.frame()
  rsds$chr <- chsq_table$chr
  rsds$celltype <- data$celltype[[1]]
  
  rsds$pvalue <- chsq_results$p.value[1]
  return(rsds)
} 

#Chisquared for each cell type 
chisq_list <- lapply(cell_type_DEG_output, csq_chr_func) %>% 
  bind_rows() %>% #collapse not_sig, SR_bias and ST_bias into single column 
  pivot_longer(c("Unbiased", "SR-biased", "ST-biased"), names_to = "Expression class", values_to = "chisq_resid")

chisq_list %>% dplyr::select(celltype, pvalue) %>% 
  unique()

#Chi squared for all cell types
dif_exp_data_table %>% group_by(`Cell type`) %>%
  summarise(`SR-biased` = sum(`SR-biased`), 
            `ST-biased` = sum(`ST-biased`),
            `Unbiased` = sum(`Unbiased`)) %>% 
  column_to_rownames("Cell type") %>%
  chisq.test(.) %>% 
  .$residuals

chisq_list %>% dplyr::select(celltype, pvalue) %>% 
  unique() %>% merge(dif_exp_data_table, by.x = 'celltype', by.y = 'Cell type') %>% 
  mutate(pvalue = format(pvalue, scientific = TRUE, digits = 3)) %>% 
  mutate(pvalue = p.adjust(pvalue)) %>%
  relocate(pvalue, .after = Unbiased) %>% 
  rename(`Cell type` = celltype, 
         `FDR Adjusted X2 p-value` = pvalue) %>% 
  .[c(3,4,7,8,5,6,1,2,9:14),] %>% 
  write.table(., "data/DEG_table.tsv", sep = "\t", quote = F, row.names = F)

##############################################################



################# INVERSION CHECK --------------
DGE_inversion <- ortholog_table[,c("REF_GENE_NAME", "start", "end")] %>% 
  merge(dif_exp_data, by.x = 'REF_GENE_NAME', by.y = 'geneID')

inversions <- data.frame(inversions = c("inv1", "inv2", "inv3", "inv4", "inv5"),
                         starts = c(76.7, 80.2, 59,32,1.2), 
                         ends = c(80.4, 87.07, 76.7, 46.3, 9.03)) %>% 
  mutate(starts = starts * 1e6, 
         ends = ends * 1e6)

inv_func <- function(x, inversions){
  if (x[12] == "Chr_X"){
    st <- inversions[,1][which(inversions[,2] < as.numeric(x[2]) & inversions[,3] > 
                                 as.numeric(x[3]))]
    if(length(st) == 0){
      return("outside inversion")
    } else {
      return(st[1])
      #  return("inversion")
    }
  }else {
    return("autosomal")
  }
}


DGE_inversion$inversion <- apply(DGE_inversion, 1, inv_func, inversions = inversions) %>% unlist()
chisq_inv_table <- DGE_inversion %>% filter(chr == "Chr_X") %>% 
  group_by(inversion, Significant, celltype) %>% 
  summarise(count = n()) %>% 
  reshape2::dcast(celltype + inversion ~ Significant, value.var = "count") 
chisq_inv_table[is.na(chisq_inv_table)] <- 0


subset_chisq <- function(cell, data){
  ss <- filter(data, celltype == cell) %>% 
    dplyr::select(-c(celltype, inversion)) %>% 
    chisq.test()
  
  return(ss)
}

inv_chisqs <- lapply(unique(chisq_inv_table$celltype), subset_chisq, data = chisq_inv_table)

#There are no interesting results for enrichment of DGE inside/outside the inversion 

###################



################### I LIKE BUT WON"T USE #####################
dif_exp_figure <- dif_exp_data %>% 
  filter(Significant != "Unbiased") %>% 
  ggplot(aes(x = celltype, fill = Significant)) + geom_bar(position = 'dodge') + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + labs(y = "Number of DEGs") 

DGE_chisq_plot <- ggplot(chisq_list, aes(x = chr, y = chisq_resid, fill = `Expression class`)) + 
  facet_wrap(~celltype) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_classic() + 
  scale_fill_brewer(palette = "Set2") + 
  labs(x = "Chromosome", y = "Residuals")

ggarrange(DGE_chisq_plot, dif_exp_figure, ncol = 1, heights = c(2, 1), common.legend = T, 
          labels = c("A", "B"), legend = "bottom")
ggsave("plots/DGE_chisq_plot.pdf", width = 15, height = 20, units = "cm")

################################################################################













########## DELETE ---------

#### CHECK CELLTYPE ABUNDANCE DIFFERENCES ----
abundances <- table(sce$celltype, sce$sample)
abundances <- unclass(abundances) 
head(abundances)

extra.info <- colData(sce)[match(colnames(abundances), sce$sample),]
y.ab <- DGEList(abundances, samples=extra.info)
keep <- filterByExpr(y.ab, group=y.ab$samples$treatment)
y.ab <- y.ab[keep,]
summary(keep)


design <- model.matrix(~factor(treatment), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)


plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)


-------------------------------------------
  
  
  
  ### OLD ABUNDANCE DIFFERENCES --------
seurat_final$broad_type <- "Somatic"
seurat_final$broad_type[seurat_final$celltype %in% c(
  "Spermatids", "Spermatocytes", "GSC/Spermatogonia")] <- "Germ"
#sce$broad_type <- sce$celltype
seurat_final_ss <- subset(seurat_final, celltype %in% c("Cyst", "GSC/Spermatogonia", "Spermatids", "Spermatocytes"))
DimPlot(seurat_final_ss, group.by = "broad_type", label = T)
cell_type_abundances_summary <- data.frame(type = seurat_final_ss$broad_type, 
                                           sample = seurat_final_ss$sample, 
                                           treatment = seurat_final_ss$treatment) %>% 
  dplyr::count(treatment, type, sample) %>% 
  spread(type,n)
cell_type_abundances_summary$prop <- cell_type_abundances_summary$Germ / 
  (cell_type_abundances_summary$Germ + cell_type_abundances_summary$Somatic)


glmer(cbind(cell_type_abundances_summary$Germ, 
            cell_type_abundances_summary$Germ + cell_type_abundances_summary$Somatic) ~ 
        cell_type_abundances_summary$treatment + (1|sample), data = cell_type_abundances_summary, family = "binomial") %>% 
  summary()


plot_data <- cell_type_abundances_summary %>% 
  mutate(Total = Germ + Somatic,
         Germ = Germ/Total, 
         Somatic = Somatic/Total) %>% 
  pivot_longer(c("Germ", "Somatic"), names_to = "Tissue",
               values_to = "Proportion of cell type")
plot_data[plot_data == 'sr'] <- "Drive"
plot_data[plot_data == 'st'] <- "Standard"
plot_data %>% 
  ggplot(aes(y = `Proportion of cell type`, x = treatment, fill = Tissue)) + geom_boxplot() + 
  theme_classic() 
plot_data %>% 
  ggplot(aes(y = `Proportion of cell type`, x = treatment, fill = Tissue)) + geom_boxplot() + 
  theme_classic() 
ggsave("plots/sex_unfolded_figures/comp_cell_number.pdf", height = 4, width = 5)


sum_data <- seurat_final_ss@meta.data %>% 
  group_by(celltype, treatment, sample) %>% 
  summarise(n = n()) %>% 
  spread(celltype,n)

sum_data$total <- rowSums(sum_data[,c(3:6)])
sum_data[,3:6] <- sum_data[,3:6]/sum_data$total

sum_data %>% 
  pivot_longer(c("Cyst", "GSC/Spermatogonia", "Spermatocytes", "Spermatids"), names_to = "celltype", values_to = "prop") %>%
  ggplot(aes(x = treatment, y = prop, fill = celltype)) + geom_boxplot() +
  theme_classic()
x <- cell_type_DEG_output$GSC_and_Spermatogonia
csq_chr_func <- function(x){
  data <- x$all
  data$sig <- 'Unbiased'
  data$sig[data$PValue < 0.05 & data$logFC > 1] <- 'ST_bias'
  data$sig[data$PValue < 0.05 & data$logFC < -1] <- 'SR_bias'
  chsq_table <- data %>%
    group_by(sig, chr) %>% 
    summarise(n = n()) %>% #trnasformn variables into rows and columns
    spread(sig, n)
  chsq_results <- chisq.test(chsq_table[,c(2:4)])
  rsds <- chsq_results$residuals %>% as.data.frame()
  rsds$chr <- chsq_table$chr
  rsds$celltype <- data$celltype[[1]]
  
  rsds$pvalue <- chsq_results$p.value[1]
  return(rsds)
} 

chisq_list <- lapply(cell_type_DEG_output, csq_chr_func) %>% 
  bind_rows() %>% #collapse not_sig, SR_bias and ST_bias into single column 
  pivot_longer(c("Unbiased", "SR_bias", "ST_bias"), names_to = "Expression class", values_to = "chisq_resid")


ggplot(chisq_list, aes(x = chr, y = chisq_resid, fill = `Expression class`)) + 
  facet_wrap(~celltype) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_classic() + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Chromosome", y = "Residuals")



