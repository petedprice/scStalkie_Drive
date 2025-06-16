#!/bin/bash
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1

cd /mnt/parscratch/users/bi1pp/Drive_sc_quant/bams/

singularity run ../sing/bcftools_v1.9-1-deb_cv1.sif bcftools mpileup -b bam.list \
	-f ../../Drive_var_rates/refs/Teleopsis_dalmanni.fna \
	-q 20 -Q 30 -C 50 -a "DP,AD" \
	| singularity run ../sing/bcftools_v1.9-1-deb_cv1.sif bcftools call -mv -Oz -o raw.vcf.gz
