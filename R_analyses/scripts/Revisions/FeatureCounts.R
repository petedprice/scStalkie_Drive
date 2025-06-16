library(edgeR)
library(tidyverse)
library(readxl)
rm(list = ls())
load("data/RData/DEG_DC.RData")
ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene, FBconcensus) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

SR_samples <- readLines('data/feature_counts/SR_bams.list') %>% 
  str_split(., pattern = "/", simplify = T) %>% .[,ncol(.)] %>% 
  gsub(".sorted.dedup.bam", "", .)
ST_samples <- readLines('data/feature_counts/ST_bams.list') %>% 
  str_split(., pattern = "/", simplify = T) %>% .[,ncol(.)] %>% 
  gsub(".sorted.dedup.bam", "", .)


write.table(c(SR_samples, ST_samples), file = 'del.csv', quote = F, row.names = F)


sample_df <- data.frame(
  Sample = c(SR_samples, ST_samples),
  Treatment = c(rep("SR", length(SR_samples)), rep("ST", length(ST_samples)))
)


fcs <- list.files("data/feature_counts/", pattern = "\\_gene_counts.txt$", full.names = T)


read_func_tmp <- function(Sample){
  file <- grep(paste0(Sample, "_gene"), fcs, value = T)
  if (length(file) ==1){
    data <- read.table(file, header = T, sep = "\t", row.names = 1)
    colnames(data)[6] <- Sample
    return(data)
  }
}

c=0
for (s in ST_samples){
  data_tmp <- read_func_tmp(s)
  data_tmp$Gene <- rownames(data_tmp)
  data_tmp <- data_tmp[,c(1,7,2,3,4,5,6)]
  if (length(data_tmp) != 0){
    if (c == 0){
      ST_data <- data_tmp
      c=c+1
    } else {
      c=c+1
      ST_data <- left_join(ST_data, data_tmp, by = c("Gene", "Chr", "Start", "End", "Strand", "Length"))
    }
  }
}

c=0
for (s in SR_samples){
  data_tmp <- read_func_tmp(s)
  data_tmp$Gene <- rownames(data_tmp)
  data_tmp <- data_tmp[,c(1,7,2,3,4,5,6)]
  
  if (length(data_tmp) != 0){
    if (c == 0){
      SR_data <- data_tmp
      c=c+1
    } else {
      SR_data <- left_join(SR_data, data_tmp, by = c("Gene", "Chr", "Start", "End", "Strand", "Length"))
    }
  }
}

colnames(SR_data)[-c(1:6)] <- paste0("SR_", colnames(SR_data)[-c(1:6)])
colnames(ST_data)[-c(1:6)] <- paste0("ST_", colnames(ST_data)[-c(1:6)])

merged <-  merge(ST_data, SR_data, by = c("Gene", "Chr", "Start", "End", "Strand", "Length"), all = T) %>% 
  `rownames<-`(.$Gene)  


counts <-merged %>%   dplyr::select(-c("Gene", "Chr", "Start", "End", "Strand", "Length"))
Jasper_CNV <- data.frame(counts = unlist(counts['STRG.3243',]), 
                         sample = c(colnames(counts))) %>% 
  mutate(treatment = str_split(sample, "_", simplify = T)[,1])

Jasper_CNV %>% 
  ggplot(aes(x = treatment, y = counts)) + geom_boxplot()

group <- str_split(colnames(counts), "_", simplify = T)[,1]

y <- DGEList(counts, genes = rownames(counts), group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- estimateDisp(y)
y <- normLibSizes(y)

fit <- glmQLFit(y, robust = T)
res <- glmQLFTest(fit)
toptags <- topTags(res)
et <- exactTest(y, pair=c("SR","ST"))

et_table <- et$table %>% 
  mutate(FDR = p.adjust(PValue, method = "fdr"), 
         Gene = rownames(et$table)) %>% 
  merge(merged[,c("Gene", "Chr")]) %>% 
  mutate(Gene = gsub("_", "-", Gene))


write.csv(et_table, 'data/CNV.csv', quote = F, row.names = F)

dif_exp_cnv <- et_table %>% 
  merge(dif_exp_data, by.x = "Gene", by.y = 'genes', 
        suffixes = c("_CNV", "_DGE"), all = T)

dif_exp_cnv <- dif_exp_cnv %>% 
  mutate(CNV_sig = 
           ifelse(FDR_CNV < 0.05 & logFC_CNV > 0, "ST", 
                  ifelse(FDR_CNV < 0.05 & logFC_CNV < 0, "SR", "UNBIAS_CNV"))) 

dif_exp_cnv %>% 
#  filter(Significant != "Unbiased") %>% 
  ggplot(aes(x = logFC_CNV, y = logFC_DGE, color = CNV_sig)) + 
  geom_point()


dif_exp_cnv %>% 
  filter(CNV_sig != "UNBIAS_CNV") %>% 
  filter(chr == "Chr_X") %>% 
  dplyr::select(Gene) %>% unique() %>% dim()



sumarised <- dif_exp_cnv %>% group_by(Gene, Chr, consensus_gene) %>% 
  summarise(ST_bias = sum(Significant == "ST-biased") > 0, 
            SR_bias = sum(Significant == "SR-biased") > 0, 
            CNV_bias = CNV_sig[1], 
            logFC_CNV = logFC_CNV[1],
            logFC_DGE = mean(logFC_DGE),
            FDR_CNV = FDR_CNV[1], 
            start = start[1]) %>% 
  mutate(DGE_bias = ifelse(ST_bias == T & SR_bias == T, "BOTH", 
                       ifelse(ST_bias == T, "ST", 
                              ifelse(SR_bias == T, "SR", "UNBIAS")))) %>% 
  dplyr::select(-c("ST_bias", "SR_bias"))
        

et_table %>% 
  filter(FDR < 0.05) %>% 
  group_by(Chr) %>% 
  summarise( n = dplyr::n())

et_table %>% 
  filter(FDR < 0.05 & abs(logFC) >= 1) %>% 
  group_by(Chr) %>% 
  summarise( n = dplyr::n())


#Filtering for genes with CNV
my_comparisons <- list(c("ST", "SR"))

#check normality 
shapiro.test(sample(sumarised$logFC_CNV, 5000, replace = F))
shapiro.test(sample(sumarised$logFC_DGE, 5000, replace = F))

A <- sumarised %>% 
  mutate(DGE_bias = factor(DGE_bias, level = c("ST", "SR", "UNBIAS"))) %>%
  #filter(CNV_bias != "UNBIAS_CNV" & DGE_bias != "UNBIAS" &  DGE_bias != "BOTH") %>% 
  filter(!is.na(DGE_bias) & !is.na(CNV_bias)) %>%
  #mutate(`DGE and coverage bias` = ifelse(CNV_bias == "SR" & DGE_bias == "SR", "SR", 
  #                       ifelse(CNV_bias == "ST" & DGE_bias == "ST", "ST", "None"))) %>% 
  #mutate(`DGE and coverage bias` = factor(`DGE and coverage bias`, level = c("ST", "SR", "None"))) %>%
  ggplot(aes(x = logFC_DGE, y= logFC_CNV)) +#, colour = `DGE and coverage bias`)) + 
  geom_point() + geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_classic() + scale_colour_brewer(palette = 'Set2') + 
  stat_cor(aes(group = 1), method = "spearman", cor.coef.name = "rho", label.x = 2, label.y = -2) + #geom_smooth(method = 'loess', se = F, colour = 'grey50', aes(group = 1)) + 
  ylab(expression(log[2]~fold~change~coverage~(ST:SR))) + 
  xlab(expression(log[2]~fold~change~DGE~(ST:SR))) + 
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5))) + 
  labs(colour = "") + 
  ylim(-3.6,2) +
  theme(legend.position = 'none')# Square legend keys 
A




B <- sumarised %>% 
  mutate(DGE_bias = factor(DGE_bias, level = c("ST", "SR", "UNBIAS"))) %>%
  filter(CNV_bias != "UNBIAS_CNV") %>% 
  filter(DGE_bias != "UNBIAS") %>% 
  filter(!is.na(DGE_bias) | !is.na(CNV_bias)) %>% 
  ggplot(aes(x = DGE_bias, y = logFC_CNV, fill = DGE_bias)) + 
  geom_boxplot() + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  labs(x = "Expression bias", fill = "") + 
 # ylab(expression(log[2]~fold~change~coverage~(ST:SR))) + 
  ylab("") +
  guides(fill = 'none') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  ylim(-3.6,2) +
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox')

  

CNV_figure <- ggarrange(A,B, labels = c("(a)", "(b)"), ncol = 2, widths = c(5,3), label.x = -0.02)

spearman_data <- sumarised %>% 
  mutate(DGE_bias = factor(DGE_bias, level = c("ST", "SR", "UNBIAS"))) %>% 
  filter(CNV_bias != "UNBIAS_CNV" & DGE_bias != "UNBIAS" &  DGE_bias != "BOTH")
shapiro.test(spearman_data$logFC_CNV)
shapiro.test(spearman_data$logFC_DGE)

dim(spearman_data)
spearman_data[abs(spearman_data$logFC_CNV) >=1,]
cor.test(spearman_data$logFC_CNV, spearman_data$logFC_DGE, method = 'spearman')

ggsave("plots/manuscript_plots/S10.tiff", height = 4, width = 6)
system("open plots/manuscript_plots/S10.tiff")


## SUP TABLE COMPARING CANDIDATES with reindhart
Rh_sup <- read_excel("data/Supplement_tables_GBE_revision_Feb14.xlsx", sheet = 5)[,c(1:10)]
colnames(Rh_sup) <- Rh_sup[1,]
Rh_sup <- Rh_sup[-1,]
Rh_ot <- merge(Rh_sup, dif_exp_cnv, by.x = "Symbol", by.y = "FBconcensus", all.x = T, all.y = F)


celltypes <- c("M", "EC", "LC", "GSC", "PS", "SS", "EST", "LST")
names(celltypes) <- c("Muscle", "Early cyst", "Late cyst", "GSC/Spermatogonia", 
                      "Primary spermatocytes", "Secondary spermatocytes", "Early spermatids", "Late spermatids")


ct_exp <- filter(dif_exp_cnv, FBconcensus %in% Rh_sup$FBID_KEY |
                   consensus_gene %in% c("STRG.3243")) %>% 
  mutate(consensus_gene = ifelse(consensus_gene == "STRG.3243", "JASPer", consensus_gene)) %>% 
  filter(!Gene %in% c("gene-7460", "gene-5012", "gene-5011")) %>% 
  mutate(combo = paste0(round(logFC_DGE,3), " (p = ", round(FDR_DGE, 3), ")")) %>% 
  dplyr::select(celltype, Gene, combo) %>% 
  pivot_wider(names_from = 'celltype', values_from = 'combo') 

colnames(ct_exp)[-1] <- paste0("logFC (fdr p-value) in ", colnames(ct_exp)[-1])
ct_exp %>% View()
sup_table_tmp1 <- filter(dif_exp_cnv, FBconcensus %in% Rh_sup$FBID_KEY |
                           consensus_gene %in% c("STRG.3243")) %>% 
  mutate(consensus_gene = ifelse(consensus_gene == "STRG.3243", "JASPer", consensus_gene)) %>% 
  filter(!Gene %in% c("gene-7460", "gene-5012", "gene-5011")) %>% 
  group_by(consensus_gene, Gene, logFC_CNV, FDR_CNV, chr, CNV_sig) %>% 
  summarise(`Expressed in` = paste(unique(celltype), collapse = ', '),
            `SR DGE in` = paste(unique(celltype[Significant == "SR-biased"]), collapse = ", "),
            `ST DGE in` = paste(unique(celltype[Significant == "ST-biased"]), collapse = ", ")) %>% 
  as.data.frame() %>% 
  dplyr::rename(Symbol = consensus_gene, 
         Ref_gene = Gene)# %>% 

sup_table_tmp1[is.na(sup_table_tmp1)] <- "Insufficient coverage"
sup_table_tmp2 <- merge(sup_table_tmp1, ct_exp, by.x = "Ref_gene", by.y = "Gene")

sup_table_tmp2 %>% 
  write.table(., file = "data/CNV_sup_table.tsv", quote = F, sep = '\t', 
              row.names = F)



### OUR TOP 6 candidates 

candidates <- filter(dif_exp_cnv, 
                     abs(logFC_CNV) >=1 & 
                       FDR_CNV < 0.05 & Significant != "Unbiased") %>% 
  dplyr::select(Gene) %>% unlist() %>% unique()
candidates


#top_tmp <- 
dif_exp_cnv %>% #mutate(DGE_bias = factor(DGE_bias, level = c("ST", "SR", "UNBIAS"))) %>%
  filter(Gene %in% candidates) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids"))) %>% 
  mutate(consensus_gene = ifelse(consensus_gene == "STRG.3243", "JASPer", consensus_gene)) %>% 
  group_by(consensus_gene, Gene, logFC_CNV, FDR_CNV, chr, CNV_sig) %>% 
  summarise(`Expressed in` = paste(unique(celltype), collapse = ', '),
            `SR DGE in` = paste(unique(celltype[Significant == "SR-biased"]), collapse = ", "),
            `ST DGE in` = paste(unique(celltype[Significant == "ST-biased"]), collapse = ", ")) %>% 
  as.data.frame() %>% 
  rename(Symbol = consensus_gene, 
         Ref_gene = Gene) %>%
  mutate(Symbol = ifelse(Symbol == Ref_gene, "No ortholog", Symbol)) %>% 
  write.table("data/MANUSCRIPT/top6_candidates.tsv", sep = '\t', quote = F, row.names = F)
  






get_ortho_group_func <- function(Gene){
  Gene2 = gsub("-", "_", Gene)
  OF_tmp1 <- OF_results[grep(Gene2, OF_results$V5),]
  genes_tmp <- str_split(OF_tmp1$V5, ",", simplify = T) %>% 
    gsub(" ", "", .)
  TF <- apply(genes_tmp, 1, function(x)return(Gene2 %in% x))
  if (sum(TF) > 0){
    
    OG <- OF_tmp1[which(TF == T),]
    dros_data <- data.frame(
      Species = "dros",
      Gene = Gene, 
      Members = str_split(OG$V4, ",", simplify = T)[1,] %>% 
        gsub(" ", "",  .), 
      Orthogroup =  OG$V2
      )
      
    stalk_data <- data.frame(
      Species = "stalk",
      Gene = Gene, 
      Members = str_split(OG$V5, ",", simplify = T)[1,] %>% 
        gsub(" ", "",  .), 
      Orthogroup =  OG$V2
    )
    outdata <- rbind(stalk_data, dros_data)
    
    outdata[outdata == ""] <- "None"
    return(outdata)
  } else {
    outdata <- data.frame(
      Gene = Gene, 
      Members = "None",
      Species = "None", 
      Orthogroup = "None"
    )
    return(outdata)
  }
}



top_tmp <- filter(dif_exp_cnv, 
                       FDR_CNV < 0.05 & Significant != "Unbiased")
  


OF_results <- read.table("data/orthology/Dmel_Tdal_Orthofinder_N0.tsv", sep = '\t')
Target_clusters <- OF_results[grep(paste(top_tmp$Gene, collapse = "|"), gsub("_", "-", OF_results$V5)),]



for (g in top_tmp$Gene){
  print(g)
  get_ortho_group_func(g)
}

outputs <- lapply(top_tmp$Gene, get_ortho_group_func) %>% 
  bind_rows()

syn_db_whole <- read.table("data/ref_files/fb_synonym_fb_2023_05.tsv", sep = '\t', fill = T, quote = "")

colnames(syn_db_whole) <- c("primary_FBid",  "organism_abbreviation",   "current_symbol",  "current_fullname",  
                            "fullname_synonym",     "symbol_synonym")

syn_db <- filter(syn_db_whole, organism_abbreviation == "Dmel") %>% #remove all rows where current_symbol ends in -RA or -RB or ]
  filter(grepl("FBgn", primary_FBid))

candidate_OGs <- outputs %>% 
  merge(syn_db, by.x = 'Members', by.y = 'primary_FBid', all.x = T, all.y = F) %>% 
  unique()

View(candidate_OGs)



#BEN THETAs
ST_theta <- read.csv("data/ben_pop_gen/ST_pergene_thetas.csv")
ST_theta <- ST_theta[,-1]
SR_theta <- read.csv("data/ben_pop_gen/SR_pergene_thetas.csv")
SR_theta <- SR_theta[,-1]

thetas <- merge(ST_theta, SR_theta, by = c("seqname", "start", "end"), suffixes = c("_ST", "_SR"))
thetas_ot <- merge(thetas, ortholog_table) %>% 
  merge(sumarised %>% dplyr::select(-c("start")), by.x = "REF_GENE_NAME", by.y = "Gene") %>% 
  mutate(taj_dif = TajimaD_ST - TajimaD_SR)

thetas_ot %>% 
  filter(DGE_bias != "BOTH") %>% 
  mutate(double_sig = ifelse(CNV_bias != "UNBIAS_CNV" & DGE_bias != "UNBIAS", "Bias", "Not")) %>% 
  ggplot(aes(x = double_sig, y = pi_SR/(end-start))) + 
  geom_boxplot()



View(thetas_ot)
