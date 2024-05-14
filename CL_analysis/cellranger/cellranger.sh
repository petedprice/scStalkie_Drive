# Alter the path for the work directory to your personal fastdata folder on the hpc, e.g /fastdata/sylvie/stalkies/wdir
#You may also need to edit the path to where you have cellranger saved on your system. I save all my software in a software folder in my home directory on the hpc 
#Command 1 will be the name of your run e.g drive1 or non_drive1
#Command 2 will be the path to your fastq folder
#Command 3 will be the path your transcriptome that you generated with mkref

#!/bin/bash

#$ -l h_rt=8:0:0

#$ -l rmem=32G

#$ -pe smp 8

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/scRNAseq/cellranger/wdir

~/software/cellranger-7.0.1/cellranger count --id=${1} \
	--fastqs=${2} \
	--transcriptome=${3}

