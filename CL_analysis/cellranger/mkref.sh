# Alter the path for the work directory to your personal fastdata folder on the hpc, e.g /fastdata/sylvie/stalkies/wdir
#You may also need to edit the path to where you have cellranger saved on your system. I save all my software in a software folder in my home directory on the hpc
# Command 1 will be a path to your reference file (fasta/fa/fna)
# Command 2 will be a path to you annotation file (gtf/gff)

#!/bin/bash

#$ -l h_rt=4:0:0

#$ -l rmem=16G

#$ -pe smp 4

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/scRNAseq/cellranger/wdir



~/software/cellranger-7.0.1/cellranger mkref \
	--genome=${3} \
	--fasta=${1} \
	--genes=${2}
