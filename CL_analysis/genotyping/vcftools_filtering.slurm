#!/bin/bash
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1

cd /mnt/parscratch/users/bi1pp/Drive_sc_quant/bams/

singularity run ../sing/vcftools_v0.1.16-1-deb_cv1.sif vcftools --gzvcf Chr_X_raw.vcf.gz \
	--maf 0.05 --max-missing 0.9 --minQ 30 --hwe 1e-6 --remove-indels \
	--recode --stdout | \
	singularity run ../sing/vcftools_v0.1.16-1-deb_cv1.sif bgzip > filtered.vcf.gz
singularity run ../sing/vcftools_v0.1.16-1-deb_cv1.sif tabix filtered.vcf.gz


