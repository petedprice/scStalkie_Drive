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
columns <- c('Symbol', 'FBID', 
             'Ave_Head', 'Ave_Female', 
             'Ave_Larva', 'Ave_Male', 'Ave_Ovaries', 
             'Ave_Testes', "Tau")

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
  exp <- exp_vals[,vals] %>% c() %>% unlist()
  x_hat <- exp/max(exp) %>% c()
  xhat_1 <- 1-x_hat %>% c()
  N=length(exp) %>% c()
  tau <- sum(xhat_1)/(N-1) %>% c()
  if (tau >= 0.9){
    ts <- tissue[which.max(x_hat)]
  } else {
    ts <- "Universal"
  }
  
  ts <- gsub("Ave_", "", ts)
  
  out <- data.frame(tau2 = tau, ts = ts)
  return(out)  
}


tau2s <- data.frame(tau2 = NA, ts = NA)
for (i in 1:nrow(tau_data_cndsd)){
  tau_tmp <- tau_function(tau_data_cndsd[i,])
  tau2s <- rbind(tau2s, tau_tmp)
}
tau2s <- tau2s[!is.na(tau2s[,1]),]
tau_data_cndsd <- cbind(tau_data_cndsd, tau2s)


ortho_XA <- ortholog_table %>% 
  dplyr::select(REF_GENE_NAME, consensus_gene, FBconcensus) %>% 
  unique() %>% 
  merge(XA, by.x = 'REF_GENE_NAME', by.y = 'genes', all.y = F) %>% 
  filter(Chr == "X")

tau_final <- ortho_XA %>% 
  merge(tau_data_cndsd, 
        by.x = 'FBconcensus', by.y = 'FBID', all.x = T) %>% 
  unique()


tau_final %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  filter(ts %in% c("Testes", "Universal")) %>% 
  ggplot(aes(x = celltype, y = XA, fill = ts, colour = ts)) + geom_boxplot() + 
  labs(x = "", y = "log(CPM)", fill = "Tau class") + 
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 7, method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
                      symnum.args = list(
                        cutpoints = c(0, 0.00001, 0.001, 0.05, 1), 
                        symbols = c("***", "**", "*", " ")), size = 7) + 
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black') + 
  theme_classic() + scale_fill_brewer(palette = 'Paired') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  theme(legend.position="top", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) +
  ylab(expression(log[2]~X:A~expression~ratio)) + 
  scale_y_continuous(breaks=c(-5,0,5)) + 
  geom_jitter()


tau_final %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  filter(ts %in% c("Testes", "Universal")) %>% 
  group_by(celltype, ts) %>% 
  summarise(mean = median(XA), 
            var = var(XA))
x1 <- tau_final %>% filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  filter(ts %in% c("Testes", "Universal")) %>% 
  filter(celltype == "Muscle" & ts == "Testes")
x2 <- tau_final %>% filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  filter(ts %in% c("Testes", "Universal")) %>% 
  filter(celltype == "Muscle" & ts == "Universal")


t.test(x1$XA, x2$XA)

#  scale_fill_manual(values ="Set2") +
  stat_compare_means(aes(label = ..p.signif..),  method = 'wilcox.test', method.args = list(alternative = 'two.sided'),
                      symnum.args = list(
                        cutpoints = c(0, 0.00001, 0.001, 0.05, 1),
                        symbols = c("***", "**", "*", " ")), size = 4)
















testis_genes <- filter(XA, celltype == "Early cyst") %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  merge(unique(ortholog_table[,c("REF_GENE_NAME", "consensus_gene")]), by.x = 'genes', by.y = "REF_GENE_NAME")



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



