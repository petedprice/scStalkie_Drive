library(tidyverse)
rm(list = ls())

load("data/RData/DEG_DC.RData")

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)
ortholog_table <- ortholog_table %>% 
  filter(chr == "Chr_X")
# read in vcf file and rename column 1 to "contig"
vcf_path="indata/PopGen/mkt_data/Chr_X_synnon_maf10.recode.vcf" # path to vcf file
vcf_file<-readLines(vcf_path)
vcf_SNPs<-data.frame(vcf_file[grep(pattern="#CHROM",vcf_file):length(vcf_file)])
vcf_SNPs <- data.frame(do.call('rbind', strsplit(as.character(vcf_SNPs[,1]),'\t',fixed=TRUE)))
colnames(vcf_SNPs) <- as.character(unlist(vcf_SNPs[1,]))
vcf_SNPs<-vcf_SNPs[-1,]
colnames(vcf_SNPs)[1]="contig"


genes <- unique(ortholog_table$REF_GENE_NAME)
vcf_SNPs_gene <- vector(length = ncol(vcf_SNPs) + 1)
for (gene in genes){
  print(gene)
  tmp <- filter(vcf_SNPs, as.numeric(POS) >= ortholog_table$start[ortholog_table$REF_GENE_NAME == gene] & 
                  as.numeric(POS) <= ortholog_table$end[ortholog_table$REF_GENE_NAME == gene])
  if (nrow(tmp) >0){
    tmp$gene <- gene
    vcf_SNPs_gene <- rbind(vcf_SNPs_gene, tmp)
  }
}

vcf_SNPs <- vcf_SNPs_gene
# read in allele frequencies (the output from vcftools is a bit unpractical, so we split the columns)
ST_freq_path="indata/PopGen/mkt_data/Chr_X_ST_freq.frq"
ST_freq<-read.delim(ST_freq_path, header=FALSE, skip=1)
colnames(ST_freq)<-c("contig","pos","n_alleles","n_chromosomes","ref","alt")

ST_freq<-separate(ST_freq,ref,c("ref_allele","ref_freq"), ":")
ST_freq<-separate(ST_freq,alt,c("alt_allele","alt_freq"), ":")
ST_freq$ST_ref_freq<-as.numeric(ST_freq$ref_freq)
ST_freq$ST_alt_freq<-as.numeric(ST_freq$alt_freq)
ST_freq$refalt <- ST_freq$ST_alt_freq * ST_freq$ST_ref_freq


SR_freq_path="indata/PopGen/mkt_data/Chr_X_SR_freq.frq"
SR_freq<-read.delim(SR_freq_path, header=FALSE, skip=1)
colnames(SR_freq)<-c("contig","pos","n_alleles","n_chromosomes","ref","alt")

SR_freq<-separate(SR_freq,ref,c("ref_allele","ref_freq"), ":") # separate columns into new columns based on delimiter ":"
SR_freq<-separate(SR_freq,alt,c("alt_allele","alt_freq"), ":")
SR_freq$SR_ref_freq<-as.numeric(SR_freq$ref_freq)
SR_freq$SR_alt_freq<-as.numeric(SR_freq$alt_freq)
SR_freq$refalt <- SR_freq$SR_alt_freq * SR_freq$SR_ref_freq

ST_freq$STpi <- ST_freq$refalt * 2
SR_freq$SRpi <- SR_freq$refalt * 2


# attach allele frequencies to SNP file
vcf_SNPs_tmp1 <- merge(vcf_SNPs, ST_freq[,c("contig", "pos", "ST_ref_freq", "ST_alt_freq", "STpi")], by.x = c("contig", "POS"), by.y = c("contig", "pos"))
vcf_SNPs_tmp2 <- merge(vcf_SNPs_tmp1, SR_freq[,c("contig", "pos", "SR_ref_freq", "SR_alt_freq", "SRpi")], by.x = c("contig", "POS"), by.y = c("contig", "pos"))

vcf_SNPs <- vcf_SNPs_tmp2
# only keep SNPs that are either synonymous or nonsynonymous.
vcf_SNPs<-vcf_SNPs[grep(pattern="synonymous_variant|missense_variant",vcf_SNPs$INFO),]
vcf_SNPs$annotation<-separate(vcf_SNPs,INFO,c("A","B","C","D","E","F","G","H","J","K","L","Z","X","N","V","Q","W","R"), "\\|")[,"B"] # random number of columns with random column names, needs to be more than the maximum amount of columns possible, you only need the second column
vcf_SNPs$contig<-factor(vcf_SNPs$contig) # reset factor levels
vcf_SNPs$gene<-factor(vcf_SNPs$gene) # reset factor levels

## count pN, pS, dN, dS

# create columns with 0s
vcf_SNPs$pN<-0
vcf_SNPs$pS<-0
vcf_SNPs$dN<-0
vcf_SNPs$dS<-0
vcf_SNPs$pN_ST<-0
vcf_SNPs$pS_ST<-0
vcf_SNPs$pN_SR<-0
vcf_SNPs$pS_SR<-0
vcf_SNPs <- filter(vcf_SNPs, !is.na(SR_ref_freq) & !is.na(SR_alt_freq) & !is.na(ST_ref_freq) & !is.na(ST_alt_freq))
                   
                   #SR_ref_freq != "NaN" | SR_alt_freq != "NaN" | ST_ref_freq != "NaN" | ST_alt_freq != "NaN")

# check the category a SNP falls into
for(i in 1:nrow(vcf_SNPs)) {
  if((vcf_SNPs$ST_ref_freq[i]==0 || vcf_SNPs$ST_alt_freq[i]==0) && (vcf_SNPs$SR_ref_freq[i]==0 || vcf_SNPs$SR_alt_freq[i]==0) && "synonymous_variant" %in% vcf_SNPs$annotation[i]) { vcf_SNPs$dS[i]<-1 }
  if((vcf_SNPs$ST_ref_freq[i]==0 || vcf_SNPs$ST_alt_freq[i]==0) && (vcf_SNPs$SR_ref_freq[i]==0 || vcf_SNPs$SR_alt_freq[i]==0) && "missense_variant" %in% vcf_SNPs$annotation[i]) { vcf_SNPs$dN[i]<-1 }
  if((vcf_SNPs$ST_ref_freq[i]!=0 && vcf_SNPs$ST_alt_freq[i]!=0) && (vcf_SNPs$SR_ref_freq[i]!=0 && vcf_SNPs$SR_alt_freq[i]!=0) && "synonymous_variant" %in% vcf_SNPs$annotation[i]) { vcf_SNPs$pS[i]<-1 }
  if((vcf_SNPs$ST_ref_freq[i]!=0 && vcf_SNPs$ST_alt_freq[i]!=0) && (vcf_SNPs$SR_ref_freq[i]!=0 && vcf_SNPs$SR_alt_freq[i]!=0) && "missense_variant" %in% vcf_SNPs$annotation[i]) { vcf_SNPs$pN[i]<-1 }
  if((vcf_SNPs$ST_ref_freq[i]==0 || vcf_SNPs$ST_alt_freq[i]==0) && (vcf_SNPs$SR_ref_freq[i]!=0 && vcf_SNPs$SR_alt_freq[i]!=0) && "synonymous_variant" %in% vcf_SNPs$annotation[i]) { vcf_SNPs$pS_ST[i]<-1 }
  if((vcf_SNPs$ST_ref_freq[i]==0 || vcf_SNPs$ST_alt_freq[i]==0) && (vcf_SNPs$SR_ref_freq[i]!=0 && vcf_SNPs$SR_alt_freq[i]!=0) && "missense_variant" %in% vcf_SNPs$annotation[i]) { vcf_SNPs$pN_ST[i]<-1 }
  if((vcf_SNPs$ST_ref_freq[i]!=0 && vcf_SNPs$ST_alt_freq[i]!=0) && (vcf_SNPs$SR_ref_freq[i]==0 || vcf_SNPs$SR_alt_freq[i]==0) && "synonymous_variant" %in% vcf_SNPs$annotation[i]) { vcf_SNPs$pS_SR[i]<-1 }
  if((vcf_SNPs$ST_ref_freq[i]!=0 && vcf_SNPs$ST_alt_freq[i]!=0) && (vcf_SNPs$SR_ref_freq[i]==0 || vcf_SNPs$SR_alt_freq[i]==0) && "missense_variant" %in% vcf_SNPs$annotation[i]) { vcf_SNPs$pN_SR[i]<-1 }
}


# make empty data.frame
MKT<-data.frame(gene=factor(),pN=numeric(), pS=numeric(), dN=numeric(), dS=numeric())

#sum pNs, pSs, dNs, and dSs for each contig across SNPs
for( i in 1:length(levels(vcf_SNPs$gene))) {
  temp<-vcf_SNPs[which(vcf_SNPs$gene==levels(vcf_SNPs$gene)[i]),]
  MKT<-rbind(MKT,data.frame(
    gene=as.character(temp$gene[1]),
    pN=sum(temp[,c("pN","pN_SR")]) ,
    pS=sum(temp[,c("pN","pS_SR")]),
    dN=sum(temp[,"dN"]),
    dS=sum(temp[,"dS"]), 
    SR_pi = mean(temp$SRpi), 
    ST_pi = mean(temp$STpi), 
    SR_piN = mean(temp$SRpi[temp$annotation == "missense_variant"], na.rm = T), 
    ST_piN = mean(temp$STpi[temp$annotation == "missense_variant"], na.rm = T),
    ST_piS = mean(temp$STpi[temp$annotation == "synonymous_variant"], na.rm = T),
    SR_piS = mean(temp$SRpi[temp$annotation == "synonymous_variant"], na.rm = T)))
}


# perform MKT test
MKT$pN.pS=MKT$pN/MKT$pS # calculate ration of nonsynonymous to synonymous polymorphisms
MKT$dN.dS=MKT$dN/MKT$dS # calculate ration of nonsynonymous to synonymous substitutions
MKT$fisher.test.P<-99  # create new column for p-values
MKT$NI <- MKT$pN.pS/MKT$dN.dS # calculate neutrality index
MKT <- MKT %>% 
  mutate(DOS = (dN/(dN + dS)) - (pN/(pN + pS)))

for(i in 1:nrow(MKT)){
  MKT$fisher.test.P[i]<-fisher.test(matrix(as.numeric(MKT[i,c(3,2,5,4)]), ncol=2))$p.value # calculate fisher exact test and copy p-value for every contig
  if((MKT$pN[i] == 0 && MKT$dN[i] == 0) || (MKT$pS[i] == 0 && MKT$dS[i] == 0) || (MKT$pS[i] == 0 && MKT$pN[i] == 0) || (MKT$dS[i] == 0 && MKT$dN[i] == 0)) { MKT$fisher.test.P[i]<-NA} # this lines assigns an NA to all p-values that are meaningless, because the contingency table was incomplete
  if(sum(as.numeric(MKT[i,c(3,2,5,4)])) < 3) { MKT$fisher.test.P[i]<-NA } # only use cases where total number of SNPs is higher than or equal to 3
}
MKT_noNAs <- MKT %>% 
  filter(!is.na(dN.dS) | is.infinite(dN.dS))
MKT_noNAs<-MKT[which(MKT$fisher.test.P != "NA"),] # remove all cases where fisher exact test is meaningless

# multiple hypothesis testing. It seems to be uncommon with MKT; probably because fisher exact produces low P-values only with much higher counts than are common for SNP data
MKT_noNAs$fisher.test.P<-p.adjust(MKT_noNAs$fisher.test.P, method = "BH") # correct p-values using Benjamini & Hochberg (1985) FDR


MKT_noNAs[sort(MKT_noNAs$dN.dS,decreasing = TRUE),] # sort by dN/dS ratio


DGE_MKT <- merge(MKT, dif_exp_data, by.x = 'gene', by.y = 'genes')

a <- DGE_MKT %>% 
  filter(celltype == "Spermatids") %>% 
  ggplot(aes(x = start, y = DOS)) + geom_point() + geom_smooth(method = "lm")#+ 
  #facet_wrap(~celltype)
b <- DGE_MKT %>% 
  filter(celltype == "Spermatids") %>% 
  ggplot(aes(x = start, y = abs(logFC), colour = Significant)) + geom_point() #+
#  facet_wrap(~celltype)
cowplot::plot_grid(a, b, ncol = 1)

DGE_MKT %>% 
  ggplot(aes(x = log(abs(logFC)), y = NI, colour = celltype)) + 
  geom_point(alpha = 0.1) + 
  facet_wrap(~celltype) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggpubr::stat_cor(method = 'pearson')


DGE_MKT %>% 
  ggplot(aes(y = DOS, x = Significant)) + geom_boxplot()

DGE_MKT %>% 
  group_by(Significant, celltype) %>% 
  summarise(mean = mean(DOS)) %>% 
  ggplot(aes(x = Significant, fill = Significant, y = mean)) + geom_bar(stat = 'identity', position = 'dodge') + 
  facet_wrap(~celltype)


DGE_MKT %>% 
  ggplot(aes(x = log(abs(logFC)), y = SR_piS)) +
  geom_point(alpha = 0.1) + facet_wrap(~celltype) + 
  geom_smooth(method = "lm") + 
  ggpubr::stat_cor(method = 'pearson')
DGE_MKT %>% 
  ggplot(aes(x = log(abs(logFC)), y = SR_pi)) +
  geom_point(alpha = 0.1) + facet_wrap(~celltype) + 
  geom_smooth(method = "lm") + 
  ggpubr::stat_cor(method = 'pearson')
DGE_MKT %>% 
  ggplot(aes(x = Significant, y = SR_piS)) + geom_boxplot() + 
  facet_wrap(~celltype)


DGE_MKT %>% 
  ggplot(aes(x = ST_pi, y = SR_pi)) + geom_point() + 
  geom_smooth(method = "lm") + 
  ggpubr::stat_cor(method = 'spearman')

vcf_SNPs %>% 
  ggplot(aes(x = POS, y = SRpi)) + 
  geom_point()


window <- 100000  
jump <- 25000
start = 1
c=0
while(start < max(vcf_SNPs$POS)){
  c=c+1
  end <- start + window -1
  temp <- vcf_SNPs %>% 
    filter(POS >= start & POS <= end) %>% 
    summarise(STpi = mean(STpi, na.rm  = T), SRpi = mean(SRpi, na.rm  = T), dN = sum(dN), 
              dS = sum(dS), pN = sum(pN), pS = sum(pS)) %>% 
    mutate(mid = start + (end - start)/2)
  start = start + jump
  if (c == 1){
    out_data <- temp
  } else {
    out_data <- rbind(out_data, temp)
  }
}
  
out_data %>% 
  ggplot(aes(x = mid, y = SRpi)) + geom_line()
