# Alter the path for the work directory to your personal fastdata folder on the hpc, e.g /fastdata/sylvie/stalkies/wdir
# Command 1 will be a path to where you want to save your reads

#!/bin/bash

#$ -l h_rt=4:0:0

#$ -l rmem=8G

#$ -pe smp 2

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/scRNAseq/cellranger/wdir

cp -r /shared/wright_lab_hpc/Shared/2022_scRNAseq_stalkies_gonad_disks/Sample_7-st1/ $1
cp -r /shared/wright_lab_hpc/Shared/2022_scRNAseq_stalkies_gonad_disks/Sample_3-sr1/ $1
