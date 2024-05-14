#!/bin/bash
# Request 5 gigabytes of real memory (mem)
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

singularity exec \
	/mnt/parscratch/users/bop20pp/Avian_scRNAseq/var_rates/sing/jupyter-scipy-notebook-latest.img \
	~/personal_git/Avian_scRNAseq/CL_analyses/nextflow/var_rates/software/OrthoFinder_source/orthofinder.py \
	-f /mnt/parscratch/users/bop20pp/Drive/ref_genome/helix/orthofinder/
	
