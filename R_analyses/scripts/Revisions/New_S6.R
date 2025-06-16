
dif_exp <- read.table("data/DEG_full_table.tsv", sep = "\t", header = T)
dif_cov <- read.table("data/CNV.csv", sep = ',', header = T) %>% 
  dplyr::select(Gene, logFC, FDR) %>%
  mutate(`Coverage bias` = ifelse(FDR < 0.05 & logFC < 0, "SR-biased", 
                ifelse(FDR < 0.05 & logFC > 0, "ST-biased", 'Unbiased'))) %>% 
  rename(`Coverage logFC` = logFC, 
         `Coverage FDR` = FDR) 


head(dif_exp)
head(dif_cov)

S6 <- dif_exp %>% 
  merge(dif_cov, by = 'Gene', all.x = T) 
S6[is.na(S6)] <- "Insufficient coverage"
colnames(S6) <- gsub("\\.", " ", colnames(S6))
S6$`Drosophila ortholog`[S6$Gene == "STRG.3243"] <- "JASPer"
write.table(S6, "data/MANUSCRIPT/S6_exp_cov.csv", sep = ',', quote = F, row.names = F)
w