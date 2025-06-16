library(tidyverse)
library(readxl)
rm(list = ls())
load("data/RData/DEG_DC.RData")

ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

Rh_sup <- read_excel("data/Supplement_tables_GBE_revision_Feb14.xlsx", sheet = 5)[,c(1:10)]
colnames(Rh_sup) <- Rh_sup[1,]
Rh_sup <- Rh_sup[-1,]

Rh_ot <- merge(Rh_sup, dif_exp_data, by.x = "FBID_KEY", by.y = "FBconcensus", all.x = T, all.y = F)

sup_table_tmp <- Rh_ot[,c("logFC", "logCPM", "FDR", "celltype", "Significant", "consensus_gene", "Symbol", "genes")]
sup_table_tmp[is.na(sup_table_tmp)] <- "Not expressed"


sup_table <- sup_table_tmp %>% 
  group_by(celltype, consensus_gene) %>%
  filter(!genes %in% c("gene-7460", "gene-5012", "gene-5011")) %>% 
  mutate(genes = ifelse(Symbol =='JASPer', 'gene-8880', genes), 
         consensus_gene = ifelse(Symbol == 'JASPer', 'JASPer', consensus_gene))
         
XA_candidates <- XA %>% group_by() %>% dplyr::select(celltype, genes, logcpm, treatment) %>% 
  filter(genes %in% sup_table$genes) %>%
  mutate(treatment = paste0("Mean_logcpm_", treatment)) %>% 
  pivot_wider(names_from = treatment, values_from = logcpm)
  
sup_table_final <- sup_table %>% 
  merge(XA_candidates, by = c("genes", "celltype"), all = T) %>% 
  dplyr::select(-c("logCPM"))
sup_table_final[is.na(sup_table_final)] <- "Not expressed"

sup_table_final %>% View()



  
XA %>% 
  filter(genes %in% sup_table$genes) %>% 
  View()

dif_exp_data <- 
  dif_exp_data %>% 
  mutate(new_sig = ifelse(logFC < 0 & FDR < 0.05, "SR-biased", 
                          ifelse(logFC > 0 & FDR < 0.05, "ST-biased", "Unbiased")))

  
dif_exp_data %>% 
  group_by(celltype, new_sig) %>% 
  summarise(n = n()) %>% #change new_sig values to columns 
  pivot_wider(names_from = new_sig, values_from = n)  %>% 
  mutate(`SR-biased` = ifelse(is.na(`SR-biased`), 0, `SR-biased`),
         `ST-biased` = ifelse(is.na(`ST-biased`), 0, `ST-biased`),
         Unbiased = ifelse(is.na(Unbiased), 0, Unbiased))


dif_exp_data %>% 
  group_by(celltype, Significant) %>% 
  summarise(n = n()) %>% #change new_sig values to columns 
  pivot_wider(names_from = Significant, values_from = n)

#%>% 
#  write.csv("data/no_FC_DGE.csv", quote = F, row.names = F)



