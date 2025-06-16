#!/bin/bash
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1

cd /mnt/parscratch/users/bi1pp/Drive_sc_quant/bams/

singularity run ../sing/vcftools_v0.1.16-1-deb_cv1.sif vcftools --gzvcf filtered.vcf.gz \
	--relatedness2 
	--out vcftools_relatedness

