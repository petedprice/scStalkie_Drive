################# load libraries and functions and data prep ###################
library(tidyverse)
library(Seurat)
rm(list = ls())

#### LOAD DATA AND PREP DATA ---
#load("data/RData/seurat_final.RData")
load("data/RData/DEG_DC.RData")

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table$REF_GENE_NAME<- gsub("gene_", "gene-", ortholog_table$REF_GENE_NAME)
ortholog_table$consensus_gene<- gsub("gene_", "gene-", ortholog_table$consensus_gene)

tau_data <- openxlsx::read.xlsx('data/Baker2016_S2_Dataset_rev.xlsx')  #
filtered_tau <- tau_data %>% 
  filter(!is.na(Symbol))%>% 
  dplyr::select(c('Symbol', 'Tau', 'FBID', 'Tissue.Specificity.I', 'Tissue.Specificity.II')) %>% 
  unique() %>% 
  group_by(Symbol) %>%
  filter(n() == 1)

table(tau_data$Tissue.Specificity.I)
table(filtered_tau$Tissue.Specificity.I)

ortho_XA <- ortholog_table %>% 
  dplyr::select(REF_GENE_NAME, consensus_gene, FBconcensus) %>% 
  unique() %>% 
  merge(XA, by.x = 'REF_GENE_NAME', by.y = 'genes', all.y = F) %>% 
  filter(Chr == "X")

tau_final <- ortho_XA %>% 
  merge(filtered_tau[,c('Symbol', 'Tau', 'FBID', 'Tissue.Specificity.I', 'Tissue.Specificity.II')], 
        by.x = 'FBconcensus', by.y = 'FBID', all.x = T) %>% 
  unique()


tau_final %>% 
  filter(Chr == "X"  & treatment == "ST") %>% 
  filter(Tissue.Specificity.I %in% c("T", "TO", "U", "L", "H")) %>% 
  ggplot(aes(x = celltype, y = XA, fill = Tissue.Specificity.I)) + geom_boxplot()#geom_jitter(alpha = 0.1) +


tau_final %>% 
  ggplot(aes(x = logcpm, y = Tau)) + geom_point() + 
  facet_wrap(~ Tissue.Specificity.II + celltype)

#####################
XA

testis_genes <- filter(XA, celltype == "Early cyst") %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  merge(unique(ortholog_table[,c("REF_GENE_NAME", "consensus_gene")]), by.x = 'genes', by.y = "REF_GENE_NAME")

XA %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  ggplot(aes(x = XA)) + geom_density() + facet_wrap(~celltype)

testis_genes %>% 
  ggplot(aes(x = logcpm, y = XA)) + geom_point()

compensated <- testis_genes %>% 
  filter(XA >= 0.2)

uncompensated <- testis_genes %>% 
  filter(XA <= -0.5) %>% 
  filter(XA >=-1.5)



ego_comp <- enrichGO(gene          = unique(compensated$consensus_gene), 
                universe = unique(testis_genes$consensus_gene),
                OrgDb         = org.Dm.eg.db, 
                keyType = 'SYMBOL',
                ont           = "ALL",
                pAdjustMethod = "bonferroni",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05, readable= TRUE)


ego_uncomp <- enrichGO(gene          = unique(uncompensated$consensus_gene), 
                     universe = unique(testis_genes$consensus_gene),
                     OrgDb         = org.Dm.eg.db, 
                     keyType = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "bonferroni",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05, readable= TRUE)


View(as.data.frame(ego_comp))
View(as.data.frame(ego_uncomp))



testis_genes <- tau_final %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST" & celltype == "Primary spermatocytes") %>% 
  unique()
  
compensated <- testis_genes %>% 
  filter(XA >= -0.2 & Tissue.Specificity.II == "T")

uncompensated <- testis_genes %>% 
  filter(XA <= -0.5 & Tissue.Specificity.II != "T")



