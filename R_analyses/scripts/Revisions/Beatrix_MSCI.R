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
library(GenomicFeatures)
library(ggrepel)
library(cowplot)
rm(list = ls())

#### LOAD DATA AND PREP DATA ---
load("data/RData/seurat_final.RData")
load("data/RData/DEG_DC.RData")
load("data/trajectory/sce_GAMed_8k.RData")
load("data/RData//traj_analyis.RData")

seurat_final <- JoinLayers(seurat_final)
Idents(seurat_final) <- "celltype"
levels(seurat_final) <- c("Muscle", "Early cyst", "Late cyst", 
                          "GSC/Spermatogonia", "Primary spermatocytes", 
                          "Secondary spermatocytes", "Early spermatids", "Late spermatids")
seurat_final@meta.data <- seurat_final@meta.data %>%  
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids")))

DefaultAssay(seurat_final) <- "RNA"

ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]


## AUTOSOMAL AND X-LINKED GENES ##
Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()

expr_matrix <- GetAssayData(seurat_final, 
                            assay = "RNA", 
                            layer = "data") %>% 
  as.matrix()# or "counts" or "scale.data"

#expr_matrix[expr_matrix==0] <- NA

A_expr <- colMeans(expr_matrix[Agenes, ], na.rm = T)
X_expr <- colMeans(expr_matrix[Xgenes, ], na.rm = T)
Dosage <- (X_expr/A_expr)
#rm(expr_matrix)
seurat_final <- AddMetaData(seurat_final, Dosage, 'Dosage')
seurat_final@meta.data <- seurat_final@meta.data %>% 
  mutate(MSCI_class = 
           ifelse(Dosage > 0.66, "CD", ifelse(
             Dosage <=0.66 & Dosage > 0.33, "ID", "R"
           )))


chisq_clusters <- c("GSC/Spermatogonia", "Primary spermatocytes", "Secondary spermatocytes", 
                    "Early spermatids", "Late spermatids")




outs <- list()
c = 0


for (clus in chisq_clusters){
  for(type in c("ID", "R")){
    c=c+1
    ss_data1 <- 
      seurat_final@meta.data %>% 
      filter((!celltype %in% chisq_clusters) | celltype == clus) %>% 
      mutate(celltype = ifelse(celltype == clus, clus, "other")) %>% 
      mutate(MSCI_class = ifelse(MSCI_class == type, type, "other")) %>% 
      group_by(MSCI_class, celltype) %>% 
      summarise(number = n()) %>% 
      group_by(celltype) %>%  
      mutate(percentage = 100*number/sum(number))
    ss_data2 <- ss_data1 %>% dplyr::select(-c(number)) %>%
      pivot_wider(names_from = c("MSCI_class"), values_from = "percentage")
    ss_data2[is.na(ss_data2)] <- 0
    
    cs <- chisq.test(ss_data2[,-1])
    out <- data.frame(clus = clus, 
                      type = type, 
                      pvalue = cs$p.value, 
                      df = cs$parameter, 
                      chisq = cs$statistic)
    outs[[c]] <- out
  }  
}

pvalues <- outs %>% 
  bind_rows() %>% 
  filter(pvalue < 0.05) %>% 
  mutate(y_position = 0.1, 
         sig = "***")
rownames(pvalues) <- NULL
pvalues

table(seurat_final@meta.data$celltype, seurat_final@meta.data$MSCI_class)



a <- seurat_final@meta.data %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  ggplot(aes(x =celltype, y = Dosage)) + 
  geom_boxplot() + 
  geom_hline(yintercept = c(1,0.66,0.33), linetype = 'dashed', colour = 'black', linewidth = 0.7) + 
  guides(fill = 'none') +
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) + 
  theme_classic() + 
  labs(x = "", y = "per cell X/A expression") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black"))

# Complete or partial dosage compensation: S0/Auto > 0.66.
# Lack of dosage compensation: S0/Auto < = 0.66 and S0/Auto > 0.33.
# Repressed: S0/Auto < = 0.33.


b <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype,  fill = MSCI_class)) + 
  geom_bar(position = 'fill') + 
  theme_classic() + 
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "", fill = "X-expression class", y = "proportion") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) +
  theme(legend.position = 'top') + 
  geom_text(data = pvalues, aes(x = clus, y = y_position, label = sig, fill = type), 
            fill = 'black', size = 5) 





plot <- ggarrange(a,b, nrow = 2, labels = c("(a)", "(b)"), label.x = -0.015)

ggsave("plots/manuscript_plots/SX_beatrix_MSCI.pdf", plot, height = 11, width = 6)
system("open plots/manuscript_plots/SX_beatrix_MSCI.pdf")



outs_df %>% 
  ggplot(aes(x = clus, fill = type, y = pvalue)) + 
  geom_point()

MSCI_chisq <- seurat_final@meta.data %>% 
  group_by(celltype) %>% 
  summarise(CD=length(which(Dosage >0.66)),
            ID=length(which(Dosage <=0.66 & Dosage >0.33)),
            R=length(which(Dosage <=0.33))) %>% 
  reframe(sum = CD+ID+R, 
          CD=CD/sum,
          ID=ID/sum,
          R=R/sum)

csq <- chisq.test(MSCI_chisq)
csq$residuals




