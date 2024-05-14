#!/bin/bash

#$ -l h_rt=2:0:0

#$ -l rmem=4G

#$ -pe smp 1

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/scRNAseq/wdir

module load apps/python/anaconda2-4.2.0

python --version

python2.7 /home/bop20pp/software/MeioticDrive2022/CL_analysis/orthologs_id/blast/tophits.py $1 $2
