#!/bin/bash

#$ -l h_rt=12:00:00

#$ -l rmem=25G

#$ -pe smp 1

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/mitohifi/wdir/

/usr/local/extras/Genomics/apps/anaconda_python/envs/mamba/bin/conda init bash

source ~/.bashrc

source activate py39

python /usr/local/extras/Genomics/bin/mitohifi.py -r $1 -f $2 -g $3 -t 1 -o 5


