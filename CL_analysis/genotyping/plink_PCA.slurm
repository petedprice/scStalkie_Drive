#!/bin/bash
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1

module load Anaconda3/2024.02-1

cd /mnt/parscratch/users/bi1pp/Drive_sc_quant/bams/

source activate plink 
plink \
	--vcf filtered.vcf.gz \
	--double-id \
	--allow-extra-chr \
	--set-missing-var-ids @:# \
	--make-bed --pca --out stalkie_plink_pca

conda deactivate
