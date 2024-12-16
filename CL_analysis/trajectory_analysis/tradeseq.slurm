#!/bin/bash
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1


ls

cd /mnt/parscratch/users/bop20pp/Drive/trajectory/
ls
singularity exec hattyoriiiiiii-rstudio-seurat_5.1.0-1.0.img Rscript /users/bop20pp/personal_git/scStalkie_Drive/R_analyses/scripts/Manuscript/Trajectory_analaysis_cluster.R $1 $2 $3
