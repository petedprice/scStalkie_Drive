#!/bin/bash

#$ -l h_rt=8:0:0

#$ -l rmem=32G

#$ -pe smp 8

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/scRNAseq/ref/orthofinder

conda activate orthofinder 
orthofinder -f $1 -d -og
