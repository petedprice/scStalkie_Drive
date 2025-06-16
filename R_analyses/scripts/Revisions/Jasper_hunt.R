Rh_sup <- readxl::read_excel("data/Supplement_tables_GBE_revision_Feb14.xlsx", sheet = 5)[,c(1:10)]
colnames(Rh_sup) <- Rh_sup[1,]
Rh_sup <- Rh_sup[-1,]


ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene, FBconcensus) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)


OF_results <- read.table("data/orthology/Dmel_Tdal_Orthofinder_N0.tsv", sep = '\t')
Target_clusters <- OF_results[grep(paste(Rh_sup$FBID_KEY, collapse = "|"), OF_results$V4),] %>% 
  merge(Rh_sup[,c("Symbol", "FBID_KEY")], by.x = 'V4', by.y = 'FBID_KEY')

jasper_cluster <- Target_clusters %>%
  filter(Symbol == "JASPer") %>% 
  dplyr::select(V5) %>% 
  str_split(",", simplify = T) %>% 
  gsub(" ", "", .) %>% 
  .[1,] %>% 
  gsub("_", "-", .)

filter(ortholog_table, REF_GENE_NAME %in% jasper_cluster)
filter(dif_exp_data, genes %in% jasper_cluster)

JASPER_TABLE <- filter(dif_exp_cnv, Gene %in% jasper_cluster) %>% 
  dplyr::select(-c('start', 'end')) %>% 
  merge(ortholog_table[,c("REF_GENE_NAME", 'start', 'end')], by.x = "Gene", by.y = "REF_GENE_NAME") %>% 
  group_by(Gene, Chr, CNV_sig, FDR_CNV, logFC_CNV, start, end) %>% 
  summarise(`Expressed in` = paste(unique(celltype), collapse = ', '),
            `SR DGE in` = paste(unique(celltype[Significant == "SR-biased"]), collapse = ", "),
            `ST DGE in` = paste(unique(celltype[Significant == "ST-biased"]), collapse = ", ")) %>% 
  rename(`Coverage bias` = `CNV_sig`, 
         `Coverage FDR` = FDR_CNV, 
         `Coverage logFC` = logFC_CNV) %>% 
  mutate(`Coverage bias` = ifelse(`Coverage bias` == "SR", "SR-biased", "Unbiased")) %>% 
  mutate(Information = ifelse(Gene == "STRG.3243", "Likely match to JASPer in Reinhardt et al. (2023)", ""))
  as.data.frame()
JASPER_TABLE[is.na(JASPER_TABLE)] <- ""
JASPER_TABLE[JASPER_TABLE == "NA"] <- ""
JASPER_TABLE

write.table(JASPER_TABLE, 'data/MANUSCRIPT/JASPER_PARALOGS.tsv', sep = '\t', quote = F, row.names = F)
