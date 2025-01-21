#!/bin/bash
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1

conda activate sac 

python /mnt/parscratch/users/bop20pp/Drive/scAlleleCount/scAlleleCount.py \
	--snps /mnt/parscratch/users/bop20pp/Drive/scAlleleCount/data/${1}_filt.recode.bed \
	--barcodes /mnt/parscratch/users/bop20pp/Drive/scAlleleCount/data/${1}_barcodes.tsv.gz \
	--output-prefix ${1} \
	--bamfile /mnt/parscratch/users/bop20pp/Drive/scAlleleCount/data/${1}_dup_NCR.bam


