#!/bin/bash
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1


ls

cd /mnt/parscratch/users/bop20pp/Drive/trajectory/
ls
singularity exec seurat_5.0.0.sif Rscript /users/bop20pp/personal_git/scStalkie_Drive/R_analyses/scripts/Manuscript/Trajectory_analaysis_cluster.R $1 $2 $3
