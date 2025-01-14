################# load libraries and functions and data prep ###################
library(tidyverse)
library(Seurat)
rm(list = ls())

#### LOAD DATA AND PREP DATA ---
#load("data/RData/seurat_final.RData")
load("data/RData/DEG_DC.RData")

#ortholog table loading and name changing
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table$REF_GENE_NAME<- gsub("gene_", "gene-", ortholog_table$REF_GENE_NAME)
ortholog_table$consensus_gene<- gsub("gene_", "gene-", ortholog_table$consensus_gene)

#tau data loading and naming
tau_data <- openxlsx::read.xlsx('data/Baker2016_S2_Dataset_rev.xlsx')  #
columns <- c('Symbol', 'FBID', 
             'Ave_Head', 'Ave_Female', 
             'Ave_Larva', 'Ave_Male', 'Ave_Ovaries', 
             'Ave_Testes', "Tau")

#condensing tau data, including averaging of RPKM genes that are duplicated
tau_data_cndsd <- tau_data[,columns] %>% 
  group_by(FBID) %>% 
  summarise(Ave_Male = mean(Ave_Male),
            Ave_Female = mean(Ave_Female),
            Ave_Larva = mean(Ave_Larva),
            Ave_Head = mean(Ave_Head),
            Ave_Ovaries = mean(Ave_Ovaries),
            Ave_Testes = mean(Ave_Testes),
            Tau = mean(Tau),
            n = n()) %>% 
  mutate(dup = ifelse(n > 1, "dup", "sing"))


#re-calc tau 
tau_function <- function(exp_vals){
  vals <- c(2:7) %>% c()
  tissue <- colnames(exp_vals)[vals] %>% c()
  tissue <- gsub("Ave_", "", tissue)
  exp <- exp_vals[,vals] %>% c() %>% unlist()
  x_hat <- exp/max(exp) %>% c()
  xhat_1 <- 1-x_hat %>% c()
  N=length(exp) %>% c()
  tau <- sum(xhat_1)/(N-1) %>% c()
  if (tau >= 0.7){
    ts <- tissue[which.max(x_hat)]
 # } else if (tau >= 0.5) {
  #  ts <- paste0("Weak ", tolower(tissue[which.max(x_hat)]))
  } else {
    ts <- "Universal"
  }
  
  
  out <- data.frame(tau2 = tau, ts = ts, testis_xhat = x_hat[6])
  return(out)  
}

# Re calculating Tau values
tau2s <- data.frame(tau2 = NA, ts = NA, testis_xhat = NA)
for (i in 1:nrow(tau_data_cndsd)){
  tau_tmp <- tau_function(tau_data_cndsd[i,])
  tau2s <- rbind(tau2s, tau_tmp)
}
tau2s <- tau2s[!is.na(tau2s[,1]),]
tau_data_cndsd <- cbind(tau_data_cndsd, tau2s)


#Mering ortholog table with XA table (XA is the dosage compensation values)
ortho_XA <- ortholog_table %>% 
  dplyr::select(REF_GENE_NAME, consensus_gene, FBconcensus) %>% 
  unique() %>% 
  merge(XA, by.x = 'REF_GENE_NAME', by.y = 'genes', all.y = F) %>% 
  filter(Chr == "X")

#Mering the above (ortho_XA) with the condensed tau data 
tau_final <- ortho_XA %>% 
  merge(tau_data_cndsd, 
        by.x = 'FBconcensus', by.y = 'FBID', all.x = T) %>% 
  unique()


tau_final %>% 
  write.csv("data/MANUSCRIPT/tau_data.csv")
