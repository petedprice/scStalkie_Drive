#!/usr/bin/env Rscript
rm(list=ls())


args = commandArgs(trailingOnly=TRUE)
vcf=args[1]
contig=args[2]
ct_length=as.numeric(args[3])
sample_info=read.csv(args[4], header = F)
snp_type=args[5]
package_dir=args[6]
gff=args[7]
fasta=args[8]


library(tidyverse)
library(PopGenome, lib = args[6])

print(vcf)
colnames(sample_info) <- c("sample", "treatment")
print(sample_info)

STs <- sample_info[sample_info$treatment == "ST",]$sample
SRs <-sample_info[sample_info$treatment == "SR", ]$sample
STs2 <- paste0(STs, ".2")
SRs2 <- paste0(SRs, ".2")

#READING IN GENOME
genome <- readVCF(vcf, numcols=10000, tid=contig, frompos=1, topos=ct_length, gffpath = gff, include.unknown=TRUE)
genom_all <- genome 
print('a')
#FIRST SLIDING WINDOW TRANSFORM


if (contig == "Chr_X"){
	populations <- list(STs, SRs)
	populations_all <- list(c(STs, SRs))
	genome <- set.populations(genome, populations, diploid = T)
	genome <- set.synnonsyn(genome, ref.chr=fasta)
	genom_all <- set.populations(genom_all, populations_all, diploid = T)
} else {
	populations <- list(STs, SRs)
	populations_all <- list(c(STs, SRs))
	genome <- set.populations(genome, populations, diploid = T)
       	genom_all <- set.populations(genom_all, populations_all, diploid = T)
}

print('read files')
## All samples pooled ##
if (ct_length > 100000){
	print('all samples pooled')
	genom_all <- sliding.window.transform(genom_all,width=100000,10000, type=2)
	print("pre-a")
	genom_all <- diversity.stats(genom_all, keep.site.info = F, pi = T)
	print("prea2")
	genom_all <- neutrality.stats(genom_all, FAST = TRUE)
	print("a")
	neut_all <- get.neutrality(genom_all)[[1]][,c("Tajima.D", "Fay.Wu.H", "Zeng.E")]#  %>% 
#  	as.data.frame()
	print("b")
	Pi_all <- as.data.frame(genom_all@Pi)
	colnames(Pi_all) <- c("Pi_all")
	Pi_all$n.sites_all <- genom_all@n.sites
	Pi_all$Pi_alls <- Pi_all$Pi_all/Pi_all$n.sites

	Pi_neut_all <- cbind(Pi_all, neut_all)
}
############################





#### gene level stats ####
print('gene level stats')
genome_genes <- splitting.data(genome, subsites	= "gene")
genome_genes <- neutrality.stats(genome_genes)
genome_genes <- MKT(genome_genes, do.fisher.test=TRUE)
MKT_values <-get.MKT(genome_genes)

ST_neutrality_stats <- get.neutrality(genome_genes)[[1]]
SR_neutrality_stats <- get.neutrality(genome_genes)[[2]]

colnames(ST_neutrality_stats) <- paste0("ST_", colnames(ST_neutrality_stats))
colnames(SR_neutrality_stats) <- paste0("SR_", colnames(SR_neutrality_stats))


all.stats5 <- cbind(ST_neutrality_stats, SR_neutrality_stats)
write.csv(all.stats5, paste0(snp_type, "_", contig, "_genelevel_PopGenome_Stats.csv"),
          quote = F, row.names = F)

write.table(MKT_values, paste0(snp_type, "_", contig, "_MKT_PopGenome_Stats.csv"))





######## ST vs SR  COMPARISONS 
if (ct_length > 100000){
	print('SR vs ST')
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
	colnames(Pi) <- c("STs_pi", "SRs_pi")
	Pi$n.sites <- genome@n.sites
	Pi$pos <- position
	Pi$STs_pis <- Pi$STs_pi/Pi$n.sites
	Pi$SRs_pis <- Pi$SRs_pi/Pi$n.sites

	all.stats2 <- cbind(all.stats1, Pi[,c(-4)], Pi_neut_all)

	#Neutrality stats
	genome <- neutrality.stats(genome, FAST=TRUE)
	ST_neutrality_stats <- get.neutrality(genome)[[1]]
	SR_neutrality_stats <- get.neutrality(genome)[[2]]

	colnames(ST_neutrality_stats) <- paste0("ST_", colnames(ST_neutrality_stats))
	colnames(SR_neutrality_stats) <- paste0("SR_", colnames(SR_neutrality_stats))

	all.stats3 <- cbind(all.stats2, ST_neutrality_stats, SR_neutrality_stats)

	#Diversity stats between
	genome <- diversity.stats.between(genome)
	dxy <- bind_cols(position, genome@nuc.diversity.between)
	colnames(dxy) <- c("Position", "dxy")

	all.stats4 <- cbind(all.stats3, dxy[,2])
	all.stats4$region <- rownames(all.stats4)

	write.csv(all.stats4, paste0(snp_type, "_", contig, "_PopGenome_Stats.csv"),
          	quote = F, row.names = F)
}
