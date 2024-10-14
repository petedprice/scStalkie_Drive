#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
vcf=args[1]
contig=args[2]
ct_length=as.numeric(args[3])
sample_info=read.csv(args[4])
snp_type=args[5]
package_dir=args[6]
gff=args[7]

library(tidyverse)
library(PopGenome, lib = args[6])


MAs <- sample_info[sample_info$Karyotype %in% c("AA") & sample_info$Sex == "M",]$run_accession
MAs1 <- paste0(MAs, "_hap1")
MAs2 <- paste0(MAs, "_hap2")

FAs <- sample_info[sample_info$Karyotype %in% c("A") & sample_info$Sex == "F",]$run_accession
FAs1 <- paste0(FAs, "_hap1")
FAs2 <- paste0(FAs, "_hap2")

MBs <- sample_info[sample_info$Karyotype %in% c("BB", "B") & sample_info$Sex == "M",]$run_accession
MBs1 <- paste0(MBs, "_hap1")
MBs2 <- paste0(MBs, "_hap2")

FBs <- sample_info[sample_info$Karyotype %in% c("BB", "B") & sample_info$Sex == "F",]$run_accession
FBs1 <- paste0(FBs, "_hap1")
FBs2 <- paste0(FBs, "_hap2")
print("doggy dogy")
MAs <-c(MAs1, MAs2)
MBs <-c(MBs1, MBs2)

Bs <- c(FBs1, FBs2,  MBs1, MBs2)[1:8]
As <- c(MAs1, MAs2, FAs1, FAs2)[1:8]

ZAs <- c(MAs1, MAs2, FAs1)[1:8]
ZBs <- c(FBs1, MBs1, MBs2)[1:8]



#READING IN GENOME
genome <- readVCF(vcf, numcols=10000, tid=contig, frompos=1, topos=ct_length, approx = TRUE)
genom_all <- genome 

#FIRST SLIDING WINDOW TRANSFORM

print("data read in")
if (contig == "NC_044241.2"){  
  print("running Z")
  print(length(ZAs))
  print(length(ZBs))
  populations <- list(ZAs, ZBs)
  populations_all <- list(c(ZAs, ZBs))
  
  outgroup <- c("LTF_hap1", "LTF_hap2")
  genome <- set.populations(genome, populations, diploid = F)
  genome <- set.outgroup(genome,new.outgroup=outgroup,diploid = F)
  genom_all <- set.outgroup(genome,new.outgroup=outgroup,diploid = F)
  genom_all <- set.populations(genom_all, populations_all, diploid = F)
  
} else {
  print("running autosomes")
  populations <- list(As, Bs)
  print(length(As))
  print(length(Bs))
  populations_all <- list(c(As, Bs))
  outgroup <- "LTF_hap1"
  genome <- set.populations(genome, populations, diploid = T)
  genome <- set.outgroup(genome,new.outgroup=outgroup,diploid = T)
  
  genom_all <- set.outgroup(genom_all,new.outgroup=outgroup,diploid = T)
  genom_all <- set.populations(genom_all, populations_all, diploid = T)
  
  
}

print("data read in")
## All samples pooled ##
genom_all <- sliding.window.transform(genom_all,width=100000,10000, type=2)
genom_all <- diversity.stats(genom_all, keep.site.info = F, pi = T)
genom_all <- neutrality.stats(genom_all, FAST = TRUE)
neut_all <- get.neutrality(genom_all)[[1]][,c("Tajima.D", "Fay.Wu.H", "Zeng.E")]  %>% 
  as.data.frame()
Pi_all <- as.data.frame(genom_all@Pi)
colnames(Pi_all) <- c("Pi_all")
Pi_all$n.sites_all <- genom_all@n.sites
Pi_all$Pi_alls <- Pi_all$Pi_all/Pi_all$n.sites

Pi_neut_all <- cbind(Pi_all, neut_all)

print('neutrality stats done')
############################

######## A AND B COMPARISONS 

genome <- sliding.window.transform(genome,width=100000,10000, type=2)

#Diversity Stats FST
genome <- F_ST.stats(genome, mode="nucleotide")
pairwise.FST <- t(genome@nuc.F_ST.pairwise)
position <- seq(from = 10000, to = nrow(pairwise.FST)*10000, by=10000)

all.stats1 <- bind_cols(position, pairwise.FST)
colnames(all.stats1) <- c("Pos", "FST")

#Diversity Stats Pi
genome <- diversity.stats(genome, keep.site.info = F, pi = T)
Pi <-  as.data.frame(genome@Pi)
colnames(Pi) <- c("Ahaps_pi", "Bhaps_pi")
Pi$n.sites <- genome@n.sites
Pi$pos <- position
Pi$Ahaps_pis <- Pi$Ahaps_pi/Pi$n.sites
Pi$Bhaps_pis <- Pi$Bhaps_pi/Pi$n.sites

all.stats2 <- cbind(all.stats1, Pi[,c(-4)], Pi_neut_all)

#Neutrality stats
genome <- neutrality.stats(genome, FAST=TRUE)
A_neutrality_stats <- get.neutrality(genome)[[1]]
B_neutrality_stats <- get.neutrality(genome)[[2]]

colnames(A_neutrality_stats) <- paste0("A_", colnames(A_neutrality_stats))
colnames(B_neutrality_stats) <- paste0("B_", colnames(B_neutrality_stats))

all.stats3 <- cbind(all.stats2, A_neutrality_stats, B_neutrality_stats)

#Diversity stats between
genome <- diversity.stats.between(genome)
dxy <- bind_cols(position, genome@nuc.diversity.between)
colnames(dxy) <- c("Position", "dxy")

all.stats4 <- cbind(all.stats3, dxy[,2])
all.stats4$region <- rownames(all.stats4)

write.csv(all.stats4, paste0("subsettedAll_PH_", snp_type, "_", contig, "_PopGenome_Stats.csv"),
          quote = F, row.names = F)
