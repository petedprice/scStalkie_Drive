#!/bin/bash
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1

module load Anaconda3/2024.02-1

source activate CNV

cd /mnt/parscratch/users/bi1pp/Drive_CNV


featureCounts -a ref_files/Teleopsis_dalmanni.gff -o feature_counts/check3_${1}_gene_counts.txt \
	-g "ID" -t "gene" \
	-p -C \
	-B \
	-T 8 \
	--primary --fraction \
	bams/${1}.sorted.dedup.bam




###	-M --fraction \
