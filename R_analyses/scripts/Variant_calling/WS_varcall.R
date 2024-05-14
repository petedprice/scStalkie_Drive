library(tidyverse)
x <- "indata/varcall/filtered_vcfs_relaxed/st1_filt.recode.vcf"
check_het_func <- function(x){
  sample <- str_split(x, "/", simplify = T) %>% c() %>% tail(.,1) %>% 
    str_split("_", simplify = T) %>% c() %>% head(.,1)
  vcf <- read.table(x) %>% as.data.frame()
  colnames(vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample)
  vcf$genotype <- str_split(vcf[,sample], ":", simplify = T) %>% .[,1]
  vcf$genotype2 <- 'hom'
  vcf$genotype2[vcf$genotype %in% c('0/1', "1/2")] <- 'het'
  vcf$genotype2[vcf$genotype %in% c('1/1', '0/0')] <- 'hom'
  vcf$genotype2[vcf$genotype == './.'] <- NA
  depths <- str_split(vcf[,sample], ":", simplify = T) %>% 
    .[,2] %>% 
    str_split(',', simplify = T) %>% 
    .[,1:2] 
  colnames(depths) <- c("refn", "altn")
  vcf <- cbind(vcf, depths)
  vcf$refn <- as.numeric(vcf$refn)
  vcf$altn <- as.numeric(vcf$altn)
  vcf_filt <- filter(vcf, is.na(genotype2) == F & (refn + altn) > 19)
  rm  <- which(vcf_filt$genotype2 == 'het' & (vcf_filt$refn < 9 |vcf_filt$altn < 9))
  vcf_filt <- vcf_filt[-rm,]
  vcf_filt$sample <- sample
  vcff_summary <- vcf_filt %>% group_by(CHROM) %>% 
    summarise(hom = length(which(genotype2 == 'hom')), 
              het = length(which(genotype2 == 'het'))) %>% 
    mutate(prop = het/(hom + het)) %>% 
    mutate(sample = sample)
  return(vcf_filt)
}
vcfs <- list.files("indata/varcall/ws_vcfs/", pattern = "recode.vcf", full.names = T)
sample_summaries <- lapply(vcfs, check_het_func) %>% 
  bind_rows()

#do a binned heterozygosity measure for each chromosome 
#so per chromosome, go along in sliding 10000bp bins and calculate the proportion of heterozygous sites
#then plot this per chromosome
#then plot the distribution of heterozygosity across the genome
sample_summaries %>% 
  group_by(CHROM, sample) %>% 
  mutate(bin = floor(POS/10000)) %>% 
  group_by(CHROM, bin, sample) %>% 
  summarise(hom = length(which(genotype2 == 'hom')), 
            het = length(which(genotype2 == 'het'))) %>% 
  mutate(prop = het/(hom + het)) %>% 
  ggplot(aes(x = bin, y = prop, colour = CHROM)) + geom_point(alpha = 0.01) + facet_wrap(~CHROM)

st1X_ss <- read.table("indata/varcall/st1_compiled.csv", 
                   sep = ",", header = T, stringsAsFactors = F) %>% 
  group_by(pos, CHROM) %>% 
  filter(CHROM == "Chr_X") %>% 
  summarise(refn_ss = sum(refn), altn_ss = sum(altn)) %>% 
  rename(pos_ss = pos) %>% 
  mutate(depth_ss = refn_ss + altn_ss, 
         het_ss = altn_ss/(refn_ss + altn_ss)) 
st1X_ws <- filter(sample_summaries, sample == 'st1' & CHROM == "Chr_X") %>% 
  dplyr::select(POS, refn, altn) %>% 
  rename(pos_ws = POS, refn_ws = refn, altn_ws = altn) %>% 
  mutate(depth_ws = refn_ws + altn_ws, 
         het_ws = altn_ws/(refn_ws + altn_ws))

st1Xcomp <- st1X_ss %>% merge(st1X_ws, by.x = "pos_ss", by.y = "pos_ws")
st1Xcomp$depth_diff <- st1Xcomp$depth_ss - st1Xcomp$depth_ws
hist(st1Xcomp$depth_diff)
head(st1Xcomp)
st1Xcomp %>% 
  ggplot(aes(x = depth_ss, y = depth_ws)) + geom_point(alpha = 0.01) 

st1Xcomp %>% 
  ggplot(aes(x = het_ss, y = het_ws)) + geom_point(alpha = 0.01)
